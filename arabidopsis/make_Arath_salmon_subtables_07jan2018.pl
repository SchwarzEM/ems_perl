#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use Scalar::Util qw(looks_like_number);

my %bigname2tag = (
    'AP1_st4_rep1_SRR064149.salmon_25oct2017'                 => 'AP1_st4_rep1',
    'AP1_st4_rep2_SRR064150.salmon_25oct2017'                 => 'AP1_st4_rep2',
    'AP1_st6.7_rep1_SRR064151.and.SRR064152.salmon_25oct2017' => 'AP1_st6.7_rep1',
    'AP1_st6.7_rep2_SRR064153.salmon_25oct2017'               => 'AP1_st6.7_rep2',
    'TRAP1_AP11.bc_3nt.exact_trim_4nt.salmon_25oct2017'       => 'TRAP1_AP11',
    'TRAP1_AP12.bc_3nt.exact_trim_4nt.salmon_25oct2017'       => 'TRAP1_AP12',
    'TRAP2_AP11.bc_3nt.exact_trim_4nt.salmon_25oct2017'       => 'TRAP2_AP11',
    'TRAP2_AP12.bc_3nt.exact_trim_4nt.salmon_25oct2017'       => 'TRAP2_AP12',
);

while (my $infile = <>) {
    chomp $infile;

    # Sample input lines:
    # /home/bioinformatics/sepal_rnaseq_dec2016/salmon/AP1_st4_rep1_SRR064149.salmon_25oct2017/quant.genes.sf
    # /home/bioinformatics/sepal_rnaseq_dec2016/salmon/TRAP1_AP11.bc_3nt.exact_trim_4nt.salmon_25oct2017/quant.genes.sf

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

