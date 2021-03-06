#!/usr/bin/env perl

# gffize_Hox_coords.pl -- Erich Schwarz <emsch@its.caltech.edu>, 6/8/2009.
# Purpose: convert Hox element coords. from supp. mat. of WBPaper00032330|pmid18981268 into a GFF2 file.

use strict;
use warnings;
use Getopt::Long;

my $element  = q{};
my $index    = 0;
my $chrom    = q{};
my $start_nt = q{};
my $end_nt   = q{};
my $suffix   = q{};
my $color    = q{};

GetOptions ( 'suffix:s' => \$suffix, 
             'color:s'  => \$color,
           );

$color = '0,0,255' if (! $color);

print 'track name=',
      '"BED track"',
      ' description=',
      '"BED format custom track example"',
      " color=$color\n",
      ;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A ((?:N|I) \d+) \s+ .+ \s ([IVX]+) : (\d+) \.\. (\d+) \s* /xms ) { 
        $element  = $1;
        $index    = 1;
        $chrom    = $2;
        $start_nt = $3;
        $end_nt   = $4;
    }
    elsif ( $input =~ /\A .+ \s ([IVX]+) : (\d+) \.\. (\d+) \s* /xms ) { 
        $chrom    = $1;
        $start_nt = $2;
        $end_nt   = $3;
        $index++;
    }
    print 'chr',      # for UCSC conventions
          "$chrom",
          "\t",
          "$start_nt",
          "\t",
          "$end_nt",
          "\t",
          "$element\.$index$suffix",
          "\n",
          ;
}

