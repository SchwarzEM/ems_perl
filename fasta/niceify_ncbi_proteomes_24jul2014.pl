#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

# Sample input header:
# >gi|390331338|ref|XP_003723254.1| PREDICTED: uncharacterized protein LOC100887917 [Strongylocentrotus purpuratus]
# >gi|571579039|ref|XP_006572568.1| PREDICTED: uncharacterized protein LOC102655962 [Apis mellifera]
# >gi|489257920|ref|WP_003165841.1| N-acetylmuramoyl-L-alanine amidase [Brevundimonas diminuta] >gi|328845812|gb|EGF95376.1| N-acetylmuramoyl-L-alanine amidase amiD [Brevundimonas diminuta ATCC 11568]
# >gi|557823168|ref|WP_023451441.1| N-acetylmuramoyl-L-alanine amidase [Asticcacaulis sp. AC460] >gi|557354055|gb|ESQ93549.1| N-acetylmuramoyl-L-alanine amidase [Asticcacaulis sp. AC460]

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A > /xms ) { 
        # Trim off NCBI's 'extra FASTA header' annotations, which are not really all that useful.
        if ( $input =~ /\A (> .+?) > /xms ) { 
            $input = $1;
        }

        # Then, insist on getting a clean result:
        if ( $input =~ /\A > (gi \| \d+ \| [a-z]+ \| ([^\|\s]+) \| [^\|\s]* \s+ .+ \[ ([A-Za-z]{4,5}) [^\[\]]* [^\s\]] \] \s*) \z/xms ) { 
            my $full_header = $1;
            my $prot        = $2;
            my $sp_init     = $3;
            # my $sp_lett3    = $4;
            my $prefix = $sp_init . q{_};
            $prot = $prefix . $prot ;
            print ">$prot  $full_header\n";
        }
        else { 
            die "Can't parse header: $input\n";
        }
    }
    else { 
        print "$input\n";
    }
}

