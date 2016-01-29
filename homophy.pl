#!/usr/bin/perl

use strict;
use warnings;
use Config::Simple;

#-------------------------------------------------------------------------------
# Config
#-------------------------------------------------------------------------------
my $conf = $ARGV[0];
my %conf = ();
my $dir  = $ARGV[1] || "output";
system("mkdir -p $dir");

if ( $conf ) {

	# check conf
	unless ( -e $conf ) {
		die "Configuration file does not exist, aborted!\n";
	}

	# parse conf
	print "Parsing configuration file $conf ...\n";
	my $para = new Config::Simple("$conf");
	%conf = $para->vars();

} else {

	die "No configuration file specified, aborted!\n";

}


#-------------------------------------------------------------------------------
# homophy
#-------------------------------------------------------------------------------

define_homology();

multiple_sequence_alignment();

phyml_tree();

clean_files();


#-------------------------------------------------------------------------------
# Subroutines
#-------------------------------------------------------------------------------
sub define_homology {

	my $query  = $conf{"query"};
	my $db_num = $conf{"db_num"};

	foreach my $i (1..$db_num) {
		my $db  = $conf{"db" . $i};
		my $out = $dir . "/" . "query-db" . $i . ".tsv";
		my $hit = $dir . "/" . "db" . $i . "_hit.fasta";
		homology_search($query, $db, $out);
		extract_seq($query, $db, $out, $hit);
		use_alias($hit, $conf{"db" . $i . "_alias"})
	}

}

sub homology_search {

	my ($query, $db, $out) = @_;
	my $add = " -e " . $conf{"e-value"} . " -F " . $conf{"filter"};

	unless ( -e $query) {
		die "Query file $query does not exist, aborted!\n";
	}

	unless ( -e $db) {
		die "Database file $db does not exist, aborted!\n";
	}

	system("formatdb -i $db -p T");

	print "Query: $query\n";
	print "Database: $db\n";
	print "Searching ...\n";
	system("blastall -p blastp -i $query -d $db -o $out -m 8" . "$add");

}

sub extract_seq {

	my ($query, $db, $out, $hit) = @_;

	print "Analyzing blast file ...\n";

	# caculate query_seq length
	my %seq_len = ();
	my $seq_id = undef;
	open (QUERY, $query) or die "Cannot open file $query: $!\n";
	while (my $line = <QUERY>) {
		chomp $line;
		if ($line =~ /^>(\S+)/) {
			$seq_id = $1;
		} else {
			$line =~ s/\s//g;
			$seq_len{$seq_id} += length($line);
		}
	}
	close QUERY;

	# extract hit_id
	my %hit_id = ();
	open (HIT, $out) or die "Cannot open file $out: $!\n";
	while (my $line = <HIT>) {
		my @w = split /\s/, $line;
		# w[0]: query id
		# w[1]: subject id
		# w[3]: alignment length
		my $min_hit = $conf{"min_hit"} * $seq_len{$w[0]};
		# remove short seq
		if ($w[3] >= $min_hit ) {
			# save nr seq id
			$hit_id{$w[1]} = 1;
		}
	}
	close HIT;

	# read db_seq
	my %db_seq = ();
	my $db_id = undef;
	open (DB, $db) or die "Cannot open file $db: $!\n";
	while (my $line = <DB>) {
		chomp $line;
		if ($line =~ /^>(\S+)/) {
			$db_id = $1;
		} else {
			$db_seq{$db_id} .= $line;
		}
	}
	close DB;

	print "Writing hit sequences ...\n";
	# write hit_seq
	open (HIT, ">$hit") or die "Cannot open file $hit: $!\n";
	foreach my $hit_id (sort keys %hit_id) {
		print HIT ">", $hit_id, "\n";
		print HIT $db_seq{$hit_id}, "\n";
	}
	close HIT;

}

sub multiple_sequence_alignment {

	my $query  = $conf{"query"};

	print "Aligning sequences ...\n";
	system("cat $query $dir/*.fas > $dir/combined.fas");
	system("clustalw2 -INFILE=$dir/combined.fas -TYPE=PROTEIN -OUTPUT=PHYLIP -OUTFILE=$dir/alignment.phy -QUIET 1>/dev/null 2>&1");

}

sub use_alias {

	my ($file, $alias) = @_;

	my %id_repo = ();
	my $i = 1;
	my $control = undef;
	open (IN, $file) or die "Cannot open file $file: $!\n";
	open (OUT, ">$dir/$alias.fas") or die "Cannot open file $dir/$alias.fas: $!\n";
	open (ALS, ">>$dir/alias.tsv") or die "Cannot open file $dir/alias.tsv: $!\n";
	while (my $line = <IN>) {
		if ($line =~ /^>(\S+)/) {
			my $full_id = $1;
			my ($id, $var) = split /\./, $full_id, 2;
			# use only one version of a protein in one locus
			unless (exists $id_repo{$id}) {
				my $alias_id = $alias . $i;
				$id_repo{$id} = $alias_id;
				print OUT ">", $alias_id, "\n";
				print ALS $alias_id, "\t", $full_id, "\n";
				$i++;
				$control = 1;
			} else {
				print ALS $id_repo{$id}, "\t", $full_id, "\n";
				$control = 0;
			}
		} else {
			print OUT $line if $control == 1;
		}
	}
	close ALS;
	close OUT;
	close IN;

}

sub phyml_tree {

	print "Running phyml ...\n";
	system("phyml -i $dir/alignment.phy -d aa -p -m WAG -v e -s NNI --no_memory_check --quiet");

}

sub clean_files {

	system("rm formatdb.log");
	system("rm $dir/query-db*.tsv");
	system("rm $dir/combined*");
	system("rm $dir/db*_hit.fasta ");

}
