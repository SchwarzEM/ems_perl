#!/usr/bin/env perl

# fix_repmodel_naming.pl -- Erich Schwarz <emsch@its.caltech.edu>, 1/29/2009.
# Purpose: fix unusable 'names' that RepeatModeler gives repeats, e.g., 'rnd-1_family-10#Unknown'.

use strict;
use warnings;
use Getopt::Long;

my $fasta  = q{};
my $prefix = q{};

GetOptions ( 'prefix:s'    => \$prefix,
             'fasta:s'     => \$fasta, );

if (! $fasta) { 
    die "Format: ./fix_repscout_naming.pl",
        " --prefix|-p [prefix, optional] --fasta|-f [FASTA]\n",
        ;
}

open my $FASTA, '<', $fasta or die "Can't open $fasta: $!";

while (my $input = <$FASTA>) { 
    if ( $input =~ /\A > ([^\s\#]+) (\#.*) \z/xms ) { 
        my $header_start  = $1;
        my $header_body   = $2;
        $input = '>' 
                 . $prefix 
                 . $header_start 
                 . '  '
                 . $header_body 
                 ;
    }
    print $input;
}

