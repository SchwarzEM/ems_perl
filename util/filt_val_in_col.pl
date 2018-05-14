#!/usr/bin/env perl

# filt_val_in_col.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/16/2010.
# Purpose: given a tab-defined column and a threshold (min. positive, or max. negative), filter output; ignores comment ('# ...') lines.

use strict;
use warnings;
use Getopt::Long;

my $column_number;
my $greater_than;
my $less_than;
my $maximum;
my $minimum;
my $help;

GetOptions ( 'column_number=i' => \$column_number,
             'greater_than:f'  => \$greater_than, 
             'less_than:f'     => \$less_than,
             'maximum|max:f'   => \$maximum,
             'minimum|min:f'   => \$minimum,
             'help'            => \$help, );

my @criteria = ($greater_than, $less_than, $maximum, $minimum);
my @defined_criteria = ();
foreach my $criterion (@criteria) {
    if (defined $criterion) { 
         push @defined_criteria, $criterion;
    }
}

if ( $help or (! defined $column_number ) or (! @defined_criteria ) ) { 
    die "\n",
        "Format: filt_val_in_col.pl\n",
        "            --column_number|-c [column of input to scan]\n",
        "            --maximum|--max   [maximum value allowed]\n",
        "            --minimum|--min   [minimum value allowed]\n",
        "            --greater_than|-g [maximum value *not* allowed]\n",
        "            --less_than|-l    [minimum value *not* allowed]\n",
        "            --help|-h\n",
        "        <input stream/files>\n",
        "\n",
        ; 
}

# Enforce human 1-based numbers, then map to a Unix-style 0-based index:
if ( $column_number < 1 ) { 
    die "Columns can number from 1 to N, not from 0, or less than 0.\n";
}
my $column_index = $column_number;
$column_index--;

while (my $input = <>) { 
    chomp $input;
    my $threshold;
    my $key_value;
    my @input_values = ();

    # To get a decision to print, the input line must pass *all* invoked criteria:
    my @print_criteria = ();

    # Test input line against all criteria:
    if ( $input !~ /\A \# /xms ) { 
        # Extract value from column:
        @input_values = split /\t/, $input;
        $key_value = $input_values[$column_index];

        # Prevent bogus scans:
        if (! ( defined $key_value ) ) {
            die "There exists no defined key value for these data in column number $column_number.\n";
        }

        # Subject the value to all invoked criteria.
        # Up to four different criteria are invokable; contradictory null-set criteria are tolerated.

        if ( defined $greater_than ) {
            if ( $key_value > $greater_than ) {
                push @print_criteria, 1;
            }
            else {
                push @print_criteria, 0;
            }
        }

        if ( defined $less_than ) {
            if ( $key_value < $less_than ) {
                push @print_criteria, 1;
            }
            else {
                push @print_criteria, 0;
            }
        }

       if ( defined $maximum ) {
            if ( $key_value <= $maximum ) {
                push @print_criteria, 1;
            }
            else {
                push @print_criteria, 0;
            }
        }

        if ( defined $minimum ) { 
            if ( $key_value >= $minimum ) { 
                push @print_criteria, 1; 
            }
            else { 
                push @print_criteria, 0;
            }
        }

        # Decide whether to print the line, based on all tests above:        
        my $print_decision = 1;
        foreach my $print_criterion (@print_criteria) {
            $print_decision = $print_decision * $print_criterion;
        }
        if ($print_decision) { 
            print "$input\n";
        }
    }
}

