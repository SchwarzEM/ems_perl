#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = "import ete3\nfrom ete3 import Tree\n\n";

while (my $input = <> ) {
    chomp $input;
    my $i          = 0;
    my $taxon_text = q{};
    if ( $input =~ /\A (\d+) \t ([^\t]+) \z/xms ) {
        $i          = $1;
        $taxon_text = $2;

        $taxon_text =~ s/\A\s+//;
        $taxon_text =~ s/\s+\z//; 

        my @taxa = map { "\"$_\"" } split '; ', $taxon_text;
        my $out_tax_text = join ', ', @taxa;

        my $out_tree = "orthofinder.$i.pruned.tree.nwk";

        print $header if $header;
        $header = q{};

        print 'tree = Tree("orthofinder.1.tree.v2.nwk")', "\n";
        print 'alignment_taxa = {', $out_tax_text, '}', "\n";
        print 'tree.prune(alignment_taxa, preserve_branch_length=True)', "\n";
        print 'tree.write(outfile="', $out_tree, '", format=1)', "\n";
        print "\n";
    }
    else {
        die "Cannot parse: $input\n";
    }
}

print "quit()\n";
