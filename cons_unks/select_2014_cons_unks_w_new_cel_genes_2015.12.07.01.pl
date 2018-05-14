#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use List::MoreUtils qw(uniq);

my %seen = ();
my $new_total = 0;

while (my $input = <>) {
    chomp $input;
    my @new_cel_genes = ();
    if ( $input =~ /c_elegans\|WBGene\d+\|\S+/xms ) {
        while ( $input =~ / (c_elegans\|WBGene\d+\|\S+) /xmsg ) {
            my $cel_gene = $1;
            $cel_gene =~ s/;\z//;
            if (! exists $seen{$cel_gene} ) {
                push @new_cel_genes, $cel_gene;
            }
        }
        if (@new_cel_genes) {
            @new_cel_genes = sort @new_cel_genes;
            @new_cel_genes = uniq @new_cel_genes;
            my $new_cel_genes = join '; ', @new_cel_genes;

            foreach my $new_cel_gene (@new_cel_genes) {
                $seen{$new_cel_gene} = 1;
            }

            $new_total    = keys %seen;
            my $new_count = @new_cel_genes;

            print "$input\t$new_cel_genes\t$new_count\t$new_total\n";
        }
    }
    else {
        if ( $input =~ /c_elegans/xms ) { 
            die "Can't parse c_elegans gene name in line: $input\n";
        }
    }
}

