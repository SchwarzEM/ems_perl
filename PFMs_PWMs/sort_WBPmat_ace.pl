#!/usr/bin/env perl

# sort_WBPmat_ace.pl -- Erich Schwarz <emsch@its.caltech.edu>, 6/8/2009.
# Purpose: given a stream of .ace matrices, print out as a single name-ordered list.

use strict;
use warnings;

my $infile     = 'infile';
my $type       = q{};
my $object     = q{};
my @temp_lines = ();

my $object2info_ref;

while (my $input = <>) { 
    if ($ARGV) { 
        $infile = $ARGV;
    }
    if ( $input =~ /\A \s* \z/xms ) { 
        @{ $object2info_ref->{$object} } = @temp_lines;
        $type       = q{};
        $object     = q{};
        @temp_lines = ();
    }

    # Order of next three 'if's matters!

    if ( ($object) and ( $input =~ /\S/xms ) ) {
        if ( $input =~ /\A Type \s+ (\S+) /xms ) {
            $type = $1;
            if ( ($type ne 'Frequency') and ($type ne 'Weight') ) { 
                die "Matrix $object is of unsupported type, $type!\n";
            }
        }
        push @temp_lines, $input;
    }
    if ( (! $object) and ( $input =~ /\S/xms ) ) { 
        if ( $input =~ /\A Position_Matrix \s+ : \s+ \" (WBPmat\d+) \"  /xms ) { 
            $object = $1;
            if ( exists $object2info_ref->{$object} ) { 
                die "Trying to enter two matrices both named $object!\n";
            }
            push @temp_lines, $input;
        }
    }
}

# Clear stored data at EOF:
@{ $object2info_ref->{$object} } = @temp_lines;
$type       = q{};
$object     = q{};
@temp_lines = ();

foreach my $matrix ( sort keys %{ $object2info_ref } ) { 
    print "\n";
    foreach my $line ( @{ $object2info_ref->{$matrix} } ) { 
        print $line;
    }
    print "\n";
}

