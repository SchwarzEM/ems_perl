#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use Scalar::Util qw(looks_like_number);

my %bigname2tag = (
    '8435_9081_52321_GCI1_may.and.june2017_ATCACG_R1.trim_filt1.salmon_08jul2017'  => 'Giant_1',
    '8435_9081_52322_GCI6_may.and.june2017_GCCAAT_R1.trim_filt1.salmon_08jul2017'  => 'Giant_6',
    '8435_9081_52323_GCI7_may.and.june2017_CAGATC_R1.trim_filt1.salmon_08jul2017'  => 'Giant_7',
    '8435_9081_52324_GCI12_may.and.june2017_CTTGTA_R1.trim_filt1.salmon_08jul2017' => 'Giant_12',
    '8435_9081_52325_SCI4_may.and.june2017_TGACCA_R1.trim_filt1.salmon_08jul2017'  => 'Small_4',
    '8435_9081_52326_SCI5_may.and.june2017_ACAGTG_R1.trim_filt1.salmon_08jul2017'  => 'Small_5',
    '8435_9081_52327_SCI8_may.and.june2017_ACTTGA_R1.trim_filt1.salmon_08jul2017'  => 'Small_8',
);

while (my $infile = <>) {
    chomp $infile;

    # Sample input line:
    # /home/bioinformatics/sepal_rnaseq_dec2016/TRAP-combined\
    # /8435_9081_52321_GCI1_may.and.june2017_ATCACG_R1.trim_filt1.salmon_08jul2017/quant.genes.sf

    if ( $infile =~ /\A \S+ \/ ([^\s\/]+) \/ quant.genes.sf \z/xms) {
        my $bigname  = $1;

        if (! exists $bigname2tag{$bigname} ) {
            die "Cannot parse big-name to tag in: $infile\n";
        }
        my $tag = $bigname2tag{$bigname};

        my $outfile  = $tag . '.quant.genes.sf.subtable';
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
                    # before rounding off decimals with 'int', add 0.5 to decimal $reads, so that values of x.5 get rounded up to x+1
                    $reads = ($reads + 0.5);
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

