#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Scalar::Util qw(looks_like_number);

my $header = "Gene\t\tAWC\tGPCR\tHouse";

while (my $input = <>) {
    chomp $input;
    # Header:
    # Gene	AWC/larvae	AWC_nz	7TM_GPCRs	Housekeeping
    # Sample data line:
    # WBGene00004487|Y57G11C.16|rps-18   2.00   7240.39   [n/a]   Housekeeping
    if ( $input =~ /\A Gene /xms ) { 
        print "$header\n" if $header;
        $header = q{};
    }
    elsif ( $input =~ /\A (WBGene\S+) \t (\S+) \t (\S+) \t ([^\t]*) \t ([^\t]*) \z/xms ) { 
        my $gene   = $1;
        my $ratio  = $2;
        my $awc_nz = $3;
        my $gpcr   = $4;
        my $hkeep  = $5;
        if ( (! looks_like_number($ratio) ) or (! looks_like_number($awc_nz) ) ) { 
            die "Can't parse putative numerical values in: $input\n";
        }
        elsif ( ( $gpcr eq '7TM_GPCRs' ) and ( $hkeep eq 'Housekeeping' ) ) { 
            # Given current data, this should not happen; but, die loudly if it does so that I can figure out how to handle it.
            die "Currently will not accept a gene that is both \"7TM_GPCRs\" and \"Housekeeping\"\n";
        }
        # Map $awc_nz to 'GPCR' data column for figure.
        elsif ( $gpcr eq '7TM_GPCRs' ) {
            print "$gene\t$ratio\t\t$awc_nz\n";
        }
        # Map $awc_nz to 'House' data column for figure.
        elsif ( $hkeep eq 'Housekeeping' ) {
            print "$gene\t$ratio\t\t\t$awc_nz\n";
        }
        else {
            # Default output line; maps $awc_nz to 'AWC' data column.
            print "$gene\t$ratio\t$awc_nz\n";
        }
    }
    else {
        die "Cannot parse input line: $input\n";
    }
}

