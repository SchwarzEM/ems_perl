#!/usr/bin/env perl

use strict;
use warnings;

my %ok_gene = ();

my $gene_list = $ARGV[0];
my $cds_file  = $ARGV[1];

open my $GENES, '<', $gene_list or die "Can't open gene list $gene_list: $!";
while (my $input = <$GENES>) { 
    chomp $input;
    if ( $input =~ /\A (Acey\S+\.t) \z/xms ) { 
        my $gene = $1;
        $ok_gene{$gene} = 1;
    } 
    else {
        die "Can't parse input line from gene list $gene_list: $input\n";
    }
}
close $GENES or die "Can't close filehandle to gene list $gene_list: $!";

open my $CDS, '<', $cds_file or die "Can't open CDS file $cds_file: $!";
while (my $input = <$CDS>) {
    chomp $input;
    if ( $input =~ /\A > ((Acey\S+\.t)\d+) /xms ) { 
        my $tx   = $1;
        my $gene = $2;
        if ($ok_gene{$gene}) { 
            print "$tx\n";
        }
    }
}
close $CDS or die "Can't close filehandle to CDS file $cds_file: $!";

