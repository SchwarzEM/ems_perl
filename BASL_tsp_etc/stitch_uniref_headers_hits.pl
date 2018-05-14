#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $uniref = $ARGV[0];
my $blastp = $ARGV[1];

my %seq2taxon = ();
my %seq2text  = ();

open my $UNIREF, '<', $uniref;
while (my $input = <$UNIREF>) {
    chomp $input;
    if ( $input =~ /\A [>] ( (\S+) .* ) \z/xms ) { 
        my $header  = $1;
        my $seqname = $2;

        if ( $header =~ /Tax=(.+) TaxID= /xms ) {
            my $taxon = $1;
            $taxon =~ s/\s+\z//;

            $seq2taxon{$seqname} = $taxon;
            $seq2text{$seqname}  = $header;
        }
        else {
            die "Cannot parse taxon ID in header: $input\n";
        }
    }
}
close $UNIREF;

open my $BLASTP, '<', $blastp;
while (my $input = <$BLASTP>) {
    chomp $input;
    if ( ( $input !~ /A [#] /xms ) and ( $input =~ /\A (\S+) \t (\S+) \t .+ \t (\S+) \t \S+ \z/xms ) ) { 
        my $query   = $1;
        my $seqhit  = $2;
        my $e_value = $3;
        if (! exists $seq2text{$seqhit} ) {
            die "Cannot find header text for \"$seqhit\" in: $input\n";
        }
        my $taxon  = $seq2taxon{$seqhit};
        my $header = $seq2text{$seqhit};
        print "$taxon\t$e_value\t$query\t$header\n";
    }
}
close $BLASTP;

