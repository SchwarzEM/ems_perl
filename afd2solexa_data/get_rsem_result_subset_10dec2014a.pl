#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use Scalar::Util qw(looks_like_number);

my $name   = $ARGV[0];
my $infile = $ARGV[1];

open my $INFILE, '<', $infile;
while (my $input = <$INFILE>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (?: \S+ \t){6} (\S+) \t \S+ \t (\S+) \t \S+ \z/xms ) { 
        my $gene  = $1;  # gene_id
        my $reads = $2;  # posterior_mean_count
        my $tpm   = $3;  # pme_TPM
        if ( $gene eq 'gene_id' ) { 
            $gene = 'Gene';
        }
        if ( $reads eq 'posterior_mean_count' ) { 
            $reads = $name . '_reads';
        }
        if ( looks_like_number($reads) ) {
            # Note that this rounds X.9 reads down to X reads; but for statistical purposes, this is acceptable, and perhaps even preferable.
            $reads = int $reads;
        }
        if ( $tpm eq 'pme_TPM' ) {
            $tpm = $name. '_TPM';
        }
        print "$gene\t$reads\t$tpm\n";
    }
    else {
        die "From input file $infile, cannot parse: $input\n";
    }
}
close $INFILE;

