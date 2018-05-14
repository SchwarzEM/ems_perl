#!/usr/bin/env perl

use strict;

my $gwas    = $ARGV[0];

open my $GWAS, '<', $gwas or die "Can't open NIH GWAS file $gwas: $!";
while (my $input = <$GWAS>) {
    chomp $input;
    if ( $input =~ /\A (?: [^\t]* \t){13} [^\t]* \t ([^\t]*) \t [^\t]* \t [^\t]* \t  /xms ) {
        my $mapped_gene = $1;
        if ( $mapped_gene =~ /\A \S+ [ ] [-] [ ] \S+ \z/xms ) { 
            my @mapped_genes = split / - /, $mapped_gene;
            foreach my $map_gene (@mapped_genes) { 
                print "$map_gene\n";
            }
        }
    }
}
close $GWAS or die "Can't close filehandle to NIH GWAS file $gwas: $!";
