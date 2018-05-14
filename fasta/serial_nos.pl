#!/usr/bin/env perl

# serial_nos.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/10/2010.
# Purpose: given a FASTA file with bad names, replace names with (optionally zero-padded) serial numbers, and displace old names to header comments.

use strict;
use warnings;
use Getopt::Long;

my $i      = 0;
my $DIGITS = 0;
my $header = q{};
my $help;

GetOptions ( 'digits=i' => \$DIGITS,
             'help'     => \$help,   );

if ($help) { 
    die "Format: serial_nos.pl -d|--digits [number of 0-padded digits in header names] -h|--help\n";
}

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A > (\S+ .*) \z/xms ) { 
        $header = $1;
        $i++;
        my $j = $i;
        if ($DIGITS) { 
            $j = sprintf("%0${DIGITS}u", $j) or die "Can't zero-pad serial number $j\n";
        }
        $input = ">$j $header";
        print "$input\n";
    }
    elsif ( $input =~ /\A > /xms ) { 
        die "Can't parse input line: $input\n";
    }
    else { 
        print "$input\n";
    }
}

