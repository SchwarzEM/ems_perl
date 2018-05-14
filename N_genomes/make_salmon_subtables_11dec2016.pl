#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use Scalar::Util qw(looks_like_number);

while (my $infile = <>) {
    chomp $infile;
    # Sample input line:
    # nigoni.vs.modENCODE.salmon.2016.11.17.01_quant.genes.sf
    if ( $infile =~ /\A \S+ \.SRR\d+ \. (\S+) _quant.genes.sf \z/xms) {
        my $tag      = $1;

        my $outfile  = $infile . '.subtable';
        $outfile     = safename($outfile);

        open my $INFILE,  '<', $infile;
        open my $OUTFILE, '>', $outfile;

        while ( my $input = <$INFILE> ) {
            chomp $input;
            # Header line:
            # Name	Length	EffectiveLength	TPM	NumReads
            if ( $input =~ /\A (\S+) (?: \t \S+){2} \t (\S+) \t (\S+) \z/xms ) {
                my $gene_id = $1;  # 'Name'
                my $tpm     = $2;  # 'TPM'
                my $reads   = $3;  # 'NumReads'

                # Change the header line so that each data column is uniquely tagged, and thus can be merged into a single large table.
                if ( $gene_id eq 'Name' ) {
                    $gene_id = 'Gene';
                    $tpm     = $tag . '_TPM'; 
                    $reads   = $tag . '_reads';
                }
                elsif ( looks_like_number($reads) ) {
                    $reads = int($reads);
                }

                # Make an output line that has only the data I want.
                my $output = "$gene_id\t$tpm\t$reads";
                print $OUTFILE "$output\n";
            }
            else { 
                die "From infile ($infile), can't parse input: $input\n";
            }
        }

        close $INFILE;
        close $OUTFILE;
    }
    else {
        die "Can't parse infile: $infile\n";
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

