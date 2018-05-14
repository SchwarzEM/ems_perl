#!/usr/bin/env perl

use strict;
use warnings;

my $targets = $ARGV[0];
my $gwas    = $ARGV[1];

my %targets = ();

my $header = q{};

open my $TARGETS, '<', $targets or die "Can't open target genes file $targets: $!";
while (my $input = <$TARGETS>) { 
    chomp $input;
    if ( $input =~ /\A (\S+) \z/xms ) { 
        my $gene = $1;
        $targets{$gene} = $1;
    }
    else { 
        die "Can't parse input line from target genes file $targets: $input\n";
    }
}
close $TARGETS or die "Can't close filehandle to target genes file $targets: $!";

open my $GWAS, '<', $gwas or die "Can't open NIH GWAS file $gwas: $!";
LOOP: while (my $input = <$GWAS>) {
    chomp $input;
    if (! $header) {
        $header = $input;
        print "$input\n";
    }
    elsif ( $input =~ /\A (?: [^\t]* \t){13} [^\t]* \t ([^\t]*) \t [^\t]* \t [^\t]* \t /xms ) {
        my $mapped_gene        = $1;
        my @mapped_genes = ();
        if ( $mapped_gene =~ /\A \S+ [ ] [-] [ ] \S+ \z/xms ) {
            @mapped_genes = split / - /, $mapped_gene;
        }
        foreach my $map_gene (@mapped_genes) { 
            if ( $targets{$map_gene} ) { 
                print "$input\n";
                next LOOP;
            }
        }
    }
    else {
        warn "Can't parse input line from NIH GWAS file $gwas: $input\n";
    }
}
close $GWAS or die "Can't close filehandle to NIH GWAS file $gwas: $!";

