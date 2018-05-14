#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use List::MoreUtils qw(uniq);

my %seen = ();
my $new_total = 0;

while (my $input = <>) {
    chomp $input;
    my @all_cel_genes          = ();
    my @new_cel_genes          = ();
    my @possible_new_cel_genes = ();

    if ( $input =~ /c_elegans\|WBGene\d+\|\S+/xms ) {
        while ( $input =~ / (c_elegans\|WBGene\d+\|\S+) /xmsg ) {
            my $cel_gene = $1;
            $cel_gene =~ s/;\z//;
            push @all_cel_genes, $cel_gene;
            # only count genes as *fully* new if they don't have CGC-style lettered names;
            # crude hack, but does save me time wasted worrying about possible ringers
            if (! exists $seen{$cel_gene} ) {
                if ( $cel_gene =~ /\A c_elegans\|WBGene\d+\|[^\|\s]+ \z/xms ) {
                    push @new_cel_genes, $cel_gene;
                }
                else {
                    push @possible_new_cel_genes, $cel_gene;
                }
            }
        }
        my $all_count = @all_cel_genes;

        if ( ( $all_count <= 2 ) and ( @new_cel_genes or @possible_new_cel_genes ) ) {
            @new_cel_genes = sort @new_cel_genes;
            @new_cel_genes = uniq @new_cel_genes;
            my $new_cel_gene_text = join '; ', @new_cel_genes;

            @possible_new_cel_genes = sort @possible_new_cel_genes;
            @possible_new_cel_genes = uniq @possible_new_cel_genes;
            my $possible_new_cel_gene_text = join '; ', @possible_new_cel_genes;


            foreach my $new_cel_gene (@new_cel_genes) {
                $seen{$new_cel_gene} = 1;
            }

            $new_total    = keys %seen;
            my $new_count = @new_cel_genes;

            print "$input\t$new_cel_gene_text\t$new_count\t$new_total\t$possible_new_cel_gene_text\n";
        }
    }
    else {
        if ( $input =~ /c_elegans/xms ) { 
            die "Can't parse c_elegans gene name in line: $input\n";
        }
    }
}

