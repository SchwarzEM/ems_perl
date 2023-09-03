#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Spec::Functions;  # catdir, catfile

my $filelist = q{};
$filelist    = $ARGV[0] if $ARGV[0];

my $i = 1;

if (! $filelist ) {
    die "Format: tab_biosample.reads_2023.09.03.01.pl [list of biosample dirs/files] > [tabulated pre-SRA data]\n";
}

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
        my $accession  = $2;
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
            print "$accession\t$condition\t$i\t$read_file1\t$read_file2\n";
            $i++;
        }
        elsif ( ( $read_file1 !~ /R1/xms ) and ( $read_file1 !~ /R2/xms ) ) {
            if (! -e $input ) {
                die "Cannot see original input file $input\n";
            }
            print "$accession\t$condition\t$i\t$read_file1\n";
            $i++;
        }
    }
    else {
        die "Cannot parse input: $input\n";
    }
}
close $FILELIST;

