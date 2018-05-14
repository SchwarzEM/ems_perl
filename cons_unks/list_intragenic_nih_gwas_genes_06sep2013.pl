#!/usr/bin/env perl

use strict;

my $gwas    = $ARGV[0];

open my $GWAS, '<', $gwas or die "Can't open NIH GWAS file $gwas: $!";
while (my $input = <$GWAS>) {
    chomp $input;
    if ( $input =~ /\A (?: [^\t]* \t){13} [^\t]* \t ([^\t]*) \t [^\t]* \t [^\t]* \t  /xms ) {
        my $mapped_gene = $1;
        if ( ( $mapped_gene !~ /\A \S+ ; \S+ \z/xms ) and ( $mapped_gene =~ /\A \S+ \z/xms ) ) {
            print "$mapped_gene\n" if ( $mapped_gene ne 'Mapped_gene' );
        }
        else {
            my @mapped_genes = ();
            if ( ( $mapped_gene =~ /\A \S+ ; \S+ \z/xms ) and ( $mapped_gene !~ /\A \S+ [ ] [-] [ ] \S+ \z/xms ) ) { 
                @mapped_genes = split /;/, $mapped_gene;
            }
            foreach my $map_gene (@mapped_genes) { 
                print "$map_gene\n";
            }
        }
    }
}
close $GWAS or die "Can't close filehandle to NIH GWAS file $gwas: $!";
