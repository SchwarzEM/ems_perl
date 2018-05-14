#!/usr/bin/env perl

# fix_repscout_naming.pl -- Erich Schwarz <emsch@its.caltech.edu>, 1/29/2009.
# Purpose: fix unusable 'names' that RepeatScout gives repeats, e.g., ">R=2 ...".

use strict;
use warnings;
use Getopt::Long;

my $fasta  = q{};
my $prefix = q{};

GetOptions ( 'prefix:s'    => \$prefix,
             'fasta:s'     => \$fasta, );

if ( (! $prefix ) and ( $fasta =~ /\w/xms ) ) { 
    $prefix = $fasta;
}

if ( (! $prefix) or (! $fasta ) ) { 
    die "Format: ./fix_repscout_naming.pl",
        " --prefix|-p [prefix, optional] --fasta|-f [FASTA]\n",
        ;
}

open my $FASTA, '<', $fasta or die "Can't open $fasta: $!";

while (my $input = <$FASTA>) { 
    if ( $input =~ /\A > (R \= (\d+)) (.*) \z/xms ) { 
        my $header_start  = $1;
        my $repeat_number = $2;
        my $header_body   = $3;
        $input = '>' 
                 . $prefix 
                 . '_Repeat_' 
                 . $repeat_number 
                 . "\t"
                 . $header_start 
                 . $header_body 
                 ;
    }
    print $input;
}

