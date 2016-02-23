#!/usr/bin/perl

use strict;
use warnings;
use Bio::DB::Taxonomy;
use Bio::Tree::Tree;
use Bio::TreeIO;

my $conf  = $ARGV[0];
my $alias_dir = $ARGV[1];

# read organisms
my $org = `grep "alias" $conf | cut -f3 -d ' ' | sed 's/\_/ /g'`;
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
		$org =~ s/ /\_/g;
		my $hit = `grep $org $alias | cut -f1`;
		$hit =~ s/$org//g;
		$org =~ s/\_/ /g;
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

# output species tree
tax2tree(@org);

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

	# based on:
	# https://github.com/gawbul/bioinformatics-scripts/blob/master/tax_identifier.pl
	foreach my $taxname (@taxnames) {
		my $taxon = $db->get_taxon(-name => $taxname);
		if (defined($taxon)) {
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
	}

	return %lineages;
}

# based on bioperl script bp_taxonomy2tree.pl
sub tax2tree {
	my @species = @_;

	my $nodesfile = "data/taxonomy/nodes.dmp";
	my $namesfile = "data/taxonomy/names.dmp";

	my $db = Bio::DB::Taxonomy->new(-source => 'flatfile',
		-nodesfile => $nodesfile,
		-namesfile => $namesfile);

	# the full lineages of the species are merged into a single tree
	my $tree = undef;
	for my $name (@species) {
		my $ncbi_id = $db->get_taxonid($name);
		if ($ncbi_id) {
			my $node = $db->get_taxon(-taxonid => $ncbi_id);
			if ($tree) {
				$tree->merge_lineage($node);
			} else {
				$tree = new Bio::Tree::Tree(-node => $node);
			}
		} else {
			warn "no NCBI Taxonomy node for species ",$name,"\n";
		}
	}

	# simple paths are contracted by removing degree one nodes
	$tree->contract_linear_paths;

	# convert tree ids to their names for nice output with TreeIO
	foreach my $node ($tree->get_nodes) {
		$node->id($node->node_name);
	}

	# the tree is output in Newick format
	my $output = new Bio::TreeIO(-format => 'newick', -file => ">species.tree");
	$output->write_tree($tree);
}
