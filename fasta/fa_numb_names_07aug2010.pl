#!/usr/bin/env perl

# fa_numb_names_07aug2010.pl -- Erich Schwarz <emsch@its.caltech.edu>, 8/7/2010.  LEGACY: kept available for backward compatibility.
# Purpose: replace FASTA names with '0000x' series; goes well with tag_FASTA_names.pl.

use strict;
use warnings;
use Getopt::Long;

my $header = q{};
my $pad    = 0;
my $i      = 1;
my $index  = q{};
my $help;

GetOptions ( 'padding=s' => \$pad,
             'help'      => \$help
           );

if ($help or ( $pad < 1) ) { 
    die "Format: fa_numb_names.pl --padding|-p [1+ digits of 0-padding] --help|-h\n";
}

while (my $input = <>) { 
    chomp $input;
    if ( $input !~ /\A > /xms ) { 
        print "$input\n";
    }
    elsif ( $input =~ /\A > (\S .*) \z/xms ) {
        $header = $1;
        my $format = '%0' . $pad . 'u';
        $index = sprintf "$format", $i;    # This lets me feed user-specified values into sprintf.
        print ">$index  $header\n";
        $i++;
    }
    else { 
        die "Can't parse input line: $input\n";
    }
}

