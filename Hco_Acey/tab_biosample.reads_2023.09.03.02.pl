#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Spec::Functions;  # catdir, catfile

my $filelist = q{};
$filelist    = $ARGV[0] if $ARGV[0];

my $header = "biosample_accession"    #   [...]
             . "\t"
             . "library_ID"           #   [replicate name]
             . "\t"
             . "title"                #   [replicate name -- can manually clarify this]
             . "\t"
             . "library_strategy"     #   RNA-Seq         # Random sequencing of whole transcriptome      
             . "\t"
             . "library_source"       #   TRANSCRIPTOMIC  # Transcription products or non genomic DNA (EST, cDNA, RT-PCR, screened libraries)
             . "\t"
             . "library_selection"    #   cDNA            # complementary DNA
             . "\t"
             . "library_layout"       #   Paired-end [or] Single
             . "\t"
             . "platform"             #   ILLUMINA
             . "\t"
             . "instrument_model"     #   unspecified     # unless detailed instrument documentation exists
             . "\t"
             . "design_description"   #   [Cite Chang et al. (2021), PubMed 32797452.]
             . "\t"
             . "filetype"             #   fastq
             . "\t"
             . "filename"             #   [file 1]
             . "\t"
             . "filename2"            #   [optional file 2]
             . "\t"
             . "filename3"            #   [ - ]
             . "\t"
             . "filename4"            #   [ - ]
             . "\t"
             . "assembly"             #   [ - ]
             . "\t"
             . "fasta_file"           #   [ - ]
             ;

if (! $filelist ) {
    die "Format: tab_biosample.reads_2023.09.03.01.pl [list of biosample dirs/files] > [tabulated pre-SRA data]\n";
}

my $biosample_accession = q{};
my $library_ID          = q{};
my $title               = q{};
my $library_strategy    = 'RNA-Seq';
my $library_source      = 'TRANSCRIPTOMIC';
my $library_selection   = 'cDNA';
my $library_layout      = q{};
my $platform            = 'ILLUMINA';
my $instrument_model    = 'unspecified';
my $design_description  = 'Single-worm mRNA-seq was performed as in Chang et al. (2021), Methods Mol. Biol. 2170, 79-99, PubMed 32797452.';
my $filetype            = 'fastq';
my $filename            = q{};
my $filename2           = q{};
my $filename3           = q{};
my $filename4           = q{};
my $assembly            = q{};
my $fasta_file          = q{};

open my $FILELIST, '<', $filelist;
while ( my $input = <$FILELIST> ) {
    chomp $input;
    if (! -r $input ) {
        die "Cannot read input file: $input\n";
    }
    # Sample inputs:
    #
    # Single read file per biological replicate:
    # SAMN37203690.WT_G15/10_WT_G15_trm.fastq.gz

    if ( $input =~ /\A ((SAMN\d+)\.([^\s\/]+)) \/ (\S+) \z/xms ) {
        my $full_name  = $1;
        my $accession = $2;
        my $condition  = $3;
        my $read_file1 = $4;
        my $read_file2 = q{};

        if ( $read_file1 =~ /R1/xms ) {
            $read_file2 = $read_file1;
            $read_file2 =~ s/R1/R2/g;

            my $orig_file = catfile($full_name, $read_file1);
            my $alt_file  = catfile($full_name, $read_file2);

            if ( $input ne $orig_file ) {
                die "Cannot reconcile original input file $input with reconstructed input file $orig_file\n";
            }
            if (! -e $input ) {
                die "Cannot see original input file $input\n";
            }
            if (! -e $alt_file ) {
                die "Cannot see second input file $alt_file\n";
            }
            if ( $input eq $alt_file ) {
                die "Cannot distinguish files $input and $alt_file\n";
            }

            # Print the header exactly once, before printing any actual data.
            print "$header\n" if $header;
            $header = q{};

            $biosample_accession = $accession;

            $library_ID = $condition;
            $title      = $condition;

            $filename  = $read_file1;
            $filename2 = $read_file2;

            $library_layout = "Paired-end";

            print "$biosample_accession\t",
                  "$library_ID\t",
                  "$title\t",
                  "$library_strategy\t",
                  "$library_source\t",
                  "$library_selection\t",
                  "$library_layout\t",
                  "$platform\t",
                  "$instrument_model\t",
                  "$design_description\t",
                  "$filetype\t",
                  "$filename\t",
                  "$filename2\t",
                  "$filename3\t",
                  "$filename4\t",
                  "$assembly\t",
                  "$fasta_file\t",
                  "\n",
                  ;
        }
        elsif ( ( $read_file1 !~ /R1/xms ) and ( $read_file1 !~ /R2/xms ) ) {
            if (! -e $input ) {
                die "Cannot see original input file $input\n";
            }

            # Print the header exactly once, before printing any actual data.
            print "$header\n" if $header;
            $header = q{};

            $biosample_accession = $accession;

            $library_ID = $condition;
            $title      = $condition;

            $filename  = $read_file1;

            $library_layout = "Single";

            print "$biosample_accession\t",
                  "$library_ID\t",
                  "$title\t",
                  "$library_strategy\t",
                  "$library_source\t",
                  "$library_selection\t",
                  "$library_layout\t",
                  "$platform\t",
                  "$instrument_model\t",
                  "$design_description\t",
                  "$filetype\t",
                  "$filename\t",
                  "$filename2\t",
                  "$filename3\t",
                  "$filename4\t",
                  "$assembly\t",
                  "$fasta_file\t",
                  "\n",
                  ;
        }
    }
    else {
        die "Cannot parse input: $input\n";
    }
}
close $FILELIST;

