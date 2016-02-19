#!/usr/bin/perl

use strict;
use warnings;

# input organism file
my $org = $ARGV[0] || "demo.organism.txt";
# match method: 0 --> exactly match
#               1 --> first word
#               2 --> first two words
#               N --> first N word(s)
my $mme = $ARGV[1] || 2;
# database directory
my $dir = $ARGV[2] || "data/seq";
# output configuration file
my $cfg = $ARGV[3] || "homophy.conf";
my $query = $ARGV[4] || "query.fa";
my $evalue = $ARGV[5] || "1e-20";
my $ncpu = $ARGV[6] || "12";

# hold organisms
my %organisms = read_organisms($org);

# hold filenames
my @filenames = ();
read_files($dir);

# make pair
my %pair = ();
foreach my $org (sort values %organisms) {
	next if (defined $pair{$org});
	foreach my $file (@filenames) {
		my @filename = split /\//, $file;
		my @orgname = split /\./, $filename[-1];
		my $name = $orgname[0];
		$name =~ s/\_/ /g;
		if ($name =~ /$org/) {
			$pair{$org} = $file;
		}
	}
	if (!defined $pair{$org}) {
		warn "$org: No matched fasta file found!\n";
	}
}

# write configuration file
my $size = keys %pair;
my $i  = 1;
open (CFG, ">$cfg") or die "Cannot open file $cfg: $!";
print CFG "# Configuration for homophy\n\n# Sequence\n";
print CFG "query  : $query\n";
print CFG "db_num : $size\n";
foreach (sort keys %pair) {
	print CFG "db" . $i . "       : " . $pair{$_} . "\n";
	s/ /\_/g;
	print CFG "db" . $i . "_alias : " . $_ . "\n";
	$i++;
}
print CFG <<BLAST;

# BLAST
e-value   : $evalue
filter    : F
min_hit   : 0.5
ncpu      : $ncpu

BLAST
close CFG;

# Get the organism name in organism file
sub read_organisms {
	my $org = shift;
	my %organisms = ();

	open (ORG, $org) or die "Cannot open file $org: $!";
	while (<ORG>) {
		chomp;
		next if /^#/;
		next if /^\S*$/;
		$organisms{$_} = $_;
		if ($mme != 0) {
			my @w = split /\s/;
			$organisms{$_} = $w[0];
			for my $id (1..($mme-1)) {
				$organisms{$_} .= " " . $w[$id];
			}
		}
		$organisms{$_} =~ s/\.//g;
	}
	close ORG;

	return %organisms;
}

# Get the full path of each file in the database directory
sub read_files {
	my $path = shift;

	opendir (DIR, $path) or die "Cannot open $path: $!";
	my @files = grep { !/^\.{1,2}$/ } readdir (DIR);
	closedir (DIR);

	@files = map { $path . '/' . $_ } @files;

	for (@files) {
		if (-d $_) {
			read_files ($_);
		} else {
			if (/\.gz$/) {
				system("gunzip -d $_");
			}
			s/\.gz$//;
			if (/\.fa$/ or /\.fas$/ or /\.fasta/ or /\.fsa/) {
				push @filenames, $_;
			}
		}
	}
}
