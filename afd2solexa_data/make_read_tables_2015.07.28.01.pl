#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use Scalar::Util qw(looks_like_number);

while (my $infile = <>) {
    chomp $infile;

    my $tag = q{};
    if ( $infile =~ /rsem_ (\S+) _2015/xms ) {
        $tag = $1;
    }
    else {
        die "Can't parse input: $infile\n";
    }

    my $outfile = 'partial_tables/' . $tag . '_reads.tsv.txt';

    $outfile = safename($outfile);

    open my $INFILE, '<', $infile;
    open my $OUTFILE, '>', $outfile;

    while (my $input = <$INFILE>) {
        chomp $input;
        if ( $input =~ /\A (\S+) (?: \s+\S+){6} \s+ (\S+) /xms ) {
            my $gene_field = $1;
            my $read_count = $2;
            my $read_header = $tag . '_reads';
            $gene_field =~ s/gene_id/Gene/;
            $read_count =~ s/posterior_mean_count/$read_header/;
            if ( looks_like_number($read_count) ) {
                $read_count = int($read_count);  # rounds partial reads down
            }
            print $OUTFILE "$gene_field\t$read_count\n";
        }
        else {
            die "In input file $infile, can't parse input line: $input\n";
        }
    }
}

sub safename {
    my $filename = $_[0];
    my $orig_filename = $filename;
    if (-e $orig_filename) {
        my $suffix1 = 1;
        $filename = $filename . ".$suffix1";
        while (-e $filename) {
            $suffix1++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix1";
        }
    }
    return $filename;
}


