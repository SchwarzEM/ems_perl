#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use Scalar::Util qw(looks_like_number);

my %orig2corr_tag = (
    ATML1__LGO_atml1_4     => 'ATML1::LGO_atml1-4',
    atml1_3                => 'atml1-3',
    lgo_2                  => 'lgo-2',
    ATML1__LGO             => 'ATML1::LGO',
    PDF1__FLAG_ATML1       => 'PDF1::FLAG-ATML1',
    PDF1__FLAG_ATML1_lgo_2 => 'PDF1::FLAG-ATML1_lgo-2',
    ATML1__LGO_atml1_3     => 'ATML1::LGO_atml1-3',
    TRAP1_AP11   => 'TRAP1_AP1.1',
    TRAP1_AP12   => 'TRAP1_AP1.2',
    TRAP1_GIANT1 => 'TRAP1_GIANT.1',
    TRAP1_GIANT2 => 'TRAP1_GIANT.2',
    TRAP1_SMALL1 => 'TRAP1_SMALL.1',
    TRAP2_AP11   => 'TRAP2_AP1.1',
    TRAP2_AP12   => 'TRAP2_AP1.2',
    TRAP2_GIANT1 => 'TRAP2_GIANT.1',
    TRAP2_GIANT2 => 'TRAP2_GIANT.2',
    TRAP2_GIANT3 => 'TRAP2_GIANT.3',
    TRAP2_SMALL1 => 'TRAP2_SMALL.1',
    TRAP2_SMALL2 => 'TRAP2_SMALL.2',
    TRAP2_SMALL3 => 'TRAP2_SMALL.3',
);

while (my $infile = <>) {
    chomp $infile;
    # Sample input line:
    # rsems_2015.09.06/rsem_atml1_3_2015.07.26.14.genes.results
    if ( $infile =~ /\A rsems_2015\.09\.06\/rsem_ (\S+) _2015\.07\.26\.\d+\.genes\.results \z/xms) {
        my $orig_tag = $1;
        my $tag      = $orig_tag;
        if ( exists $orig2corr_tag{$orig_tag} ) {
            $tag = $orig2corr_tag{$orig_tag};
        }

        my $outfile  = $infile;
        $outfile     =~ s/\.genes\.results/.genes.tpm.txt/;
        $outfile     =~ s/rsems_2015.09.06/rsems_2015.09.06_subtables/;
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
                my $output = "$gene_id\t$tpm";
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

