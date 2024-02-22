#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $infile = <>) {
    chomp $infile;
    # Sample input: /home/ems/wallacei/salmon/c_briggsae.WS275.XO.1_salmon_2020.04.01/quant.genes.sf
    # Also want to be able to process inputs with: ../quant.sf

    if ( $infile =~ / ([^\s\/]+) (?:_|.) salmon (?: (?:_|.) \d{4}\.\d{2}\.\d{2}\.\d{2,3}){0,1} \/ quant (?: \.genes){0,1} \.sf \z/xms ) {
        my $stem   = $1;

        my $outfile = "$stem.salmon.tsv.txt";
        $outfile    = safename($outfile);

        my $header = "Gene\t$stem.TPM\t$stem.reads\n";

        open my $INFILE,  '<', $infile;
        open my $OUTFILE, '>', $outfile;
        while (my $input = <$INFILE> ) {
            chomp $input;
            if ( $input =~ /\A (\S+) \t [^\t]* \t [^\t]* \t (\S+) \t (\S+) \z/xms ) {
                my $gene  = $1;
                my $tpm   = $2;
                my $reads = $3;

                if ( $gene ne 'Name' ) {
                    # round off decimal fractions in reads *after* adding 0.5, so that fractions from 0.5001 to 0.9999 round up.
                    # do not try to round 'NumReads' header text, though!
                    $reads = ($reads + 0.5);
                    $reads = int($reads);

                    print $OUTFILE $header if $header;
                    $header = q{};
                    print $OUTFILE "$gene\t$tpm\t$reads\n";
                }
            }
            else {
                die "From input file $infile, cannot parse input line: $input\n";
            }
        }
        close $INFILE;
        close $OUTFILE;
    }
    else {
        die "Cannot parse input: $infile\n";
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

