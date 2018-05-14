#!/usr/bin/env perl

# serialize_redundant_names_v01.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/10/2010.
# Purpose: given a FASTA file with redundant names, append (possibly zero-padded) serial numbers.

use strict;
use warnings;
use Getopt::Long;

my $i      = 0;
my $DIGITS = 0;
my $header = q{};
my $name   = q{};
my %name2no = ();
my $help;

GetOptions ( 'digits=i' => \$DIGITS,
             'help'     => \$help,   );

if ($help) { 
    die "Format: serial_nos.pl -d|--digits [number of 0-padded digits in header names] -h|--help\n";
}

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A > ( (\S+) .*) \z/xms ) { 
        $header = $1;
        $name   = $2;

        # First numerical value is always '1'
        my $i = $name2no{$name};
        $i++;

        # Each succeeding value is '2', '3', etc.
        $name2no{$name} = $i;

        # After incrementing number, optionally reformat how it's printed:
        if ($DIGITS) { 
            $i = sprintf("%0${DIGITS}u", $i) or die "Can't zero-pad serial number $i\n";
        }

        $input = '>' . $name . '_' . "$i  $header";
        print "$input\n";
    }
    elsif ( $input =~ /\A > /xms ) { 
        die "Can't parse input line: $input\n";
    }
    else { 
        print "$input\n";
    }
}

