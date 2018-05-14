#!/usr/bin/env perl

# bedify_gff_coords.pl -- Erich Schwarz <emsch@its.caltech.edu>, 6/9/2009.
# Purpose: extract coordinates from a GFF; put into semblance of bed format

use strict;
use warnings;
use Getopt::Long;

my $chr      = q{};
my $start_nt = q{};
my $end_nt   = q{};
my $id       = q{};
my $color    = q{};

GetOptions ( 'color:s'  => \$color, );

$color = '0,0,255' if (! $color);

print 'track name=',
      '"BED track"',
      ' description=',
      '"BED format custom track example"',
      " color=$color\n",
      ;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ / \A
                     (\S+)                    # DNA
                     \t [^\t]* \t [^\t]* \t
                     (\d+)                    # start nt
                     \t
                     (\d+)                    # end nt
                     \t [^\t]* \t [^\t]* \t [^\t]* \t
                     \S+ \s+ 
                     (\S+)                    # name (probably double-quoted)
                     .*
                     /xms ) { 
        $chr      = $1;
        $start_nt = $2;
        $end_nt   = $3;
        $id       = $4;
        $id =~ s/\"//g;
        if ( $chr !~ /\A chr/xms ) { 
            $chr = 'chr' . $chr;
        }
        print "$chr\t$start_nt\t$end_nt\t$id\n";
    }
}

