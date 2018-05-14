#!/usr/bin/env perl

use strict;
use warnings;

my %genes = ();

while (my $input = <>) { 
    chomp $input;
    $input =~ s/\AP3_//;
    if ( $input =~ /\A >/xms ) { 
        if ( $input =~ /\A > (\S+) \.t \d+ \z/xms ) { 
            my $id = $1;
            $genes{$id}++;
        }
        if ( $input !~ /\A > \S+ \.t \d+ \z/xms ) { 
            die "Couldn't parse: $input\n";
        }
    }
}

my @gene_list = sort keys %genes;
my $gene_count = scalar(@gene_list);
print "$gene_count genes\n";

my $i = 1;
my $j = 0;
while ($j < $gene_count ) { 
    my @genes_w_i = grep { $genes{$_} == $i } @gene_list;
    my $num_genes_w_i = scalar(@genes_w_i);
    print "$num_genes_w_i genes predicted to have $i isoforms\n";
    $i++;
    $j += $num_genes_w_i;
}

print "Done.\n";

