#!/usr/bin/env perl

# pick_seqnames_by_len.pl -- Erich Schwarz <emsch@its.caltech.edu>, 8/22/2008.
# Purpose: given --max=L or --min=L, export sorted names from seqnames_by_length.pl output.

use strict;
use warnings;
use Getopt::Long;

my $MINIMUM = 0;
my $MAXIMUM = 0;
my %stored  = ();

GetOptions ( 
             "lowest=i"  => \$MINIMUM, 
             "highest=i" => \$MAXIMUM, ); 

if ( ( $MINIMUM != 0 ) and ( $MAXIMUM != 0 ) and ( $MINIMUM > $MAXIMUM ) ) { 
    die "Minimum length $MINIMUM is greater than maximum length $MAXIMUM!\n";
}
if ( (! $MINIMUM != 0 ) and (! $MAXIMUM != 0 ) ) { 
    die 'Format:',
        ' ./pick_faseqs_by_len.pl',
        '  -l|--low X [and/or]  -h|--high Y',
        '  [seqnames_by_length.pl output]',
        "\n";
}

while (my $input = <>) { 
    chomp $input;
    if ($input =~ /\A (\S+) \t (\d+) \z/xms) { 
        my $name = $1;
        my $len  = $2;

        # Default to acceptance ...
        $stored{$name} = 1;

        # ... but then filter for min, max, or both:
        if ( ( $MINIMUM != 0 ) and ( $len < $MINIMUM ) ) { 
            delete $stored{$name};
        }
        if ( ( $MAXIMUM != 0 ) and ( $len > $MAXIMUM ) ) {
            delete $stored{$name};
        }
    }
    if ( $input !~ /\A \S+ \t \d+ \z/xms ) { 
        warn "Malformed input line: $input\n";
    }
}

foreach my $name (sort keys %stored) { 
    print "$name\n";
}

