#!/usr/bin/perl

use strict;
use warnings;
use Bio::DB::Taxonomy;
use Bio::Tree::Tree;

my $conf  = $ARGV[0];
my $alias_dir = $ARGV[1];

# read organisms
my $org = `grep "alias" $conf | cut -f3 -d ' ' | cut -f1,2 -d "_"`;
my @org = split /\n/, $org;

# read alias
my @alias = ();
read_alias($alias_dir);
@alias = sort by_strnum @alias;

# header
my $header = "";
foreach my $alias (@alias) {
	my @a = split /\//, $alias;
	$header .= "\t$a[-2]";
}
print $header, "\tTax\n";

# tax info
my %tax = taxonomy(@org);

# stat
my %num = ();
foreach my $org (@org) {
	foreach my $alias (@alias) {
		my $hit = `grep $org $alias | cut -f1`;
		$hit =~ s/$org//g;
		if ($hit ne "") {
			my @v = split /\n/, $hit;
			@v = sort by_num @v;
			$num{$org} .= "\t" . $v[-1];
		} else {
			$num{$org} .= "\t" . 0;
		}
	}
	print $org, $num{$org}, "\t", $tax{$org}, "\n";
}

# sort by number
sub by_num {
	$a <=> $b;
}

sub by_strnum {
	$a =~ /(\d+)/;
	my $numa = $1;
	$b =~ /(\d+)/;
	my $numb = $1;
	return $numa <=> $numb;
}
# read alias in the subdirectories
sub read_alias {
	my $path = shift;

	opendir (DIR, $path) or die "Cannot open $path: $!";
	my @files = grep { !/^\.{1,2}$/ } readdir (DIR);
	closedir (DIR);

    @files = map { $path . '/' . $_ } @files;

	for (@files) {
		if (-d $_) {
			read_alias($_);
		} else {
			if (/alias\.tsv$/) {
				push @alias, $_;
			}
		}
	}
}

# extract taxonomy information
sub taxonomy {
	my @taxnames = @_;
	my %lineages = ();

	my $nodesfile = "data/taxonomy/nodes.dmp";
	my $namesfile = "data/taxonomy/names.dmp";

	my $db = Bio::DB::Taxonomy->new(-source => 'flatfile',
		-nodesfile => $nodesfile,
		-namesfile => $namesfile);

	foreach my $taxname (@taxnames) {
		my $taxon = $db->get_taxon(-name => $taxname);
		my $tree = Bio::Tree::Tree->new(-node => $taxon);
		my @taxa = $tree->get_nodes;
		my $nodes = "";
		for (my $i = 1; $i < @taxa - 1; $i++) {
			my $t = $taxa[$i];
			my $nname = $db->get_taxon(-taxonid => $t->id());
			$nodes .= $nname->scientific_name() . "|";
		}
		$lineages{$taxname} = $nodes . $taxon->scientific_name();
	}

	return %lineages;
}
