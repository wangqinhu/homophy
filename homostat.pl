#!/usr/bin/perl

use strict;
use warnings;

my $conf  = $ARGV[0];
my $alias_dir = $ARGV[1];

# read organisms
my $org = `grep "alias" $conf | cut -f3 -d ' ' | cut -f1,2 -d "_"`;
my @org = split /\n/, $org;

# read alias
my @alias = ();
read_alias($alias_dir);

# header
my $header = "";
foreach my $alias (@alias) {
	my @a = split /\//, $alias;
	$header .= "\t$a[-2]";
}
print $header, "\n";

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
	print $org, $num{$org}, "\n";
}

# sort by number
sub by_num {
	$a <=> $b;
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
