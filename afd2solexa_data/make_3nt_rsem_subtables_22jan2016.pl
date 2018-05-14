#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use Scalar::Util qw(looks_like_number);

my %orig2corr_tag = (
    atml1_3_rep1                => 'atml1-3_rep1',
    atml1_3_rep2                => 'atml1-3_rep2',
    atml1_3_rep3                => 'atml1-3_rep3',
    atml1_4_rep1                => 'atml1-4_rep1',
    ATML1__LGO_atml1_3_rep1     => 'ATML1::LGO_atml1-3_rep1',
    ATML1__LGO_atml1_3_rep2     => 'ATML1::LGO_atml1-3_rep2',
    ATML1__LGO_atml1_3_rep3     => 'ATML1::LGO_atml1-3_rep3',
    ATML1__LGO_atml1_4_rep1     => 'ATML1::LGO_atml1-4_rep1',
    ATML1__LGO_rep1             => 'ATML1::LGO_rep1',
    ATML1__LGO_rep2             => 'ATML1::LGO_rep2',
    ATML1__LGO_rep3             => 'ATML1::LGO_rep3',
    Col_WT_rep1                 => 'Col_WT_rep1',
    Col_WT_rep2                 => 'Col_WT_rep2',
    Col_WT_rep3                 => 'Col_WT_rep3',
    lgo_2_rep1                  => 'lgo-2_rep1',
    lgo_2_rep2                  => 'lgo-2_rep2',
    lgo_2_rep3                  => 'lgo-2_rep3',
    PDF1__FLAG_ATML1_lgo_2_rep1 => 'PDF1::FLAG-ATML1_lgo-2_rep1',
    PDF1__FLAG_ATML1_lgo_2_rep2 => 'PDF1::FLAG-ATML1_lgo-2_rep2',
    PDF1__FLAG_ATML1_lgo_2_rep3 => 'PDF1::FLAG-ATML1_lgo-2_rep3',
    PDF1__FLAG_ATML1_rep1       => 'PDF1::FLAG-ATML1_rep1',
    PDF1__FLAG_ATML1_rep2       => 'PDF1::FLAG-ATML1_rep2',
    PDF1__FLAG_ATML1_rep3       => 'PDF1::FLAG-ATML1_rep3',
);

while (my $infile = <>) {
    chomp $infile;
    # Sample input line:
    # RSEM_3nt_results_2016.01.21/rsem_atml1_3_rep1.trim_exact_3nt_2016.01.21.01.genes.results
    if ( $infile =~ /\A RSEM_3nt_results_2016\.01\.21\/rsem_ (\S+) \.trim_exact_3nt_2016\.01\.21\.\d{2}\.genes\.results \z/xms) {
        my $orig_tag = $1;
        my $tag      = $orig_tag;

        if (! exists $orig2corr_tag{$orig_tag} ) {
            die "Cannot translate original tag \"$orig_tag\", and thus cannot parse input file: $infile\n";
        }
        if ( exists $orig2corr_tag{$orig_tag} ) {
            $tag = $orig2corr_tag{$orig_tag};
        }

        my $outfile  = $infile . '.subtable';
        $outfile     =~ s/RSEM_3nt_results_2016.01.21/RSEM_3nt_results_2016.01.21_subtables/;
        $outfile     = safename($outfile);

        open my $INFILE,  '<', $infile;
        open my $OUTFILE, '>', $outfile;

        while ( my $input = <$INFILE> ) {
            chomp $input;
            if ( $input =~ /\A (\S+) (?: \t \S+){6} \t (\S+) (?: \t \S+) \t (\S+) (?: \t \S+) \t (\S+) (?: \t \S+){3} \z/xms ) {
                my $gene_id = $1;  # gene_id
                my $reads   = $2;  # posterior_mean_count
                my $tpm     = $3;  # pme_TPM
                my $min_tpm = $4;  # TPM_ci_lower_bound

                # Change the header line so that each data column is uniquely tagged, and thus can be merged into a single large table.
                if ( $gene_id eq 'gene_id' ) {
                    $gene_id = 'Gene';
                    $reads   = $tag . '_reads';
                    $tpm     = $tag . '_TPM';
                    $min_tpm = $tag . '_minTPM';
                }
                elsif ( looks_like_number($reads) ) {
                    $reads = int($reads);
                }

                # Make an output line that has only the data I want.
                my $output = "$gene_id\t$reads\t$tpm\t$min_tpm";
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

