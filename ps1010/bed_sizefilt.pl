#!/usr/bin/env perl 

# bed_sizefilt.pl -- Erich Schwarz <emsch@its.caltech.edu>, 12/9/2009.
# Purpose: select minimum (and/or maximum) sizes of elements in a BED file.

use strict;
use warnings;
use Getopt::Long;

my $min;
my $max;
my $help;

GetOptions ( 'min:i' => \$min,
             'max:i' => \$max, 
             'help'  => \$help, );

if ( $help or ( (! defined $min) and (! defined $max) ) ) { 
    die "Format: bed_sizefilt.pl",
        " --min [X] --max [Y]",
        " (X, Y integers of nt; at least one is required)",
        "\n",
        ;
}

# Error-check minimum and maximum (if any were provided):
if ( ( defined $min ) and ( ( $min < 0 ) or ( $min != int($min) ) ) ) { 
    die "Minimum value $min should be a nonnegative integer!\n";
}
if ( ( defined $max ) and ( ( $max < 1 ) or ( $max != int($max) ) ) ) {
    die "Maximum value $max should be a positive integer!\n";
}
if ( ( defined $min ) and ( defined $max ) and ($min > $max) ) { 
    die "Minimum value $min should be less than or equal to",
        " maximum value $max!\n",
        ;
}

# Allow at most one header line free printing.
my $print_header = 1;

while (my $input = <>) { 
    chomp $input;

    # Print the first line if it *is* a plausible header:
    if ( $print_header and ( $input !~ /\A chr\S+ \t (\d+) \t (\d+) \s* /xms ) ) { 
        print "$input\n";
    }
    # Then, shut off that option.
    $print_header = 0;

    # For rest of input, require size-selection:
    if ( $input =~ /\A chr\S+ \t (\d+) \t (\d+) \s* /xms ) { 
        my ($nt1, $nt2, $length);
        $nt1 = $1;
        $nt2 = $2;
        if ( $nt1 > $nt2 ) { 
            die "Reversed order of nt coordinates in: $input\n";
        }
        $length = $nt2 - $nt1;  # N.B.: the actual site starts at $nt1 + 1.
        if ( ( (! defined $min ) or ( ( defined $min ) and ( $length >= $min ) ) )
             and
             ( (! defined $max ) or ( ( defined $max ) and ( $length <= $max ) ) ) ) {
            print "$input\n";
        }
    }
}

