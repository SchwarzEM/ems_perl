#!/usr/bin/env perl

# pick_minmax_column_val.pl -- Erich Schwarz <emsch@its.caltech.edu>, 12/6/2010.
# Purpose: from a tab-defined input file or stream, given N columns, leave all *other* columns intact; but, for the numbered columns, pick and print only the minimum or maximum value per line, at the end of the line; pass on comment ('# ...') llines untouched.

use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);
use List::Util qw(min max);

my $infile               = q{};

my $minimum;
my $maximum;

my $plain_value;
my $abs_value;
my $ratio;

my @input_column_numbers = ();
my @human_column_numbers = ();
my @perl_column_numbers  = ();
my $help;

GetOptions ( 'column_numbers=s{,}' => \@input_column_numbers,
             'infile=s'            => \$infile,
             'plain_value'         => \$plain_value,
             'abs_value'           => \$abs_value,
             'ratio'               => \$ratio,
             'minimum|min'         => \$minimum,
             'maximum|max'         => \$maximum,
             'help'                => \$help, );

# Default to previous behavior:
if ( (! $plain_value ) and (! $ratio ) and (! $abs_value ) ) { 
    $plain_value = 1;
}

if (    $help 
     or (! $infile                      ) 
     or (! @input_column_numbers        ) 

     or ( $plain_value and $abs_value )
     or ( $plain_value and $ratio )
     or ( $abs_value   and $ratio )

     or ( (! $minimum) and (! $maximum) ) 
     or ( $minimum and $maximum         ) ) { 
    die "Format: pick_minmax_column_val.pl\n",
        "    --infile|-i <input stream/files>\n",
        "    --minimum|--min <or> --maximum|--max  [pick one criterion]\n",
        "    --plain_value|-p      [look for largest/smallest overt value; default]\n",
        "      <or> --abs_value|-a [look for largest/smallest abs. value, with 0.0 minimum]\n",
        "      <or> --ratio|-r     [look for largest/smallest ratio, with 1.0 minimum]\n",
        "    --column_numbers|-c [columns of input to choose between]\n",
        "    --help|-h\n",
        ;
}

# Accept either a stream from '-' or a standard file.
my $INPUT_FILE;
if ($infile eq '-') {
    # Special case: get the stdin handle
    $INPUT_FILE = *STDIN{IO};
}
else {
    # Standard case: open the file
    open $INPUT_FILE, '<', $infile or die "Can't open input file $infile. $!\n";
}

# Interpret the *strings* given to --columan_numbers, so that ranges like '1-4' or '10-7' can be interpreted.
foreach my $column_number ( @input_column_numbers ) {
    my $val1 = q{};
    my $val2 = q{};
    my @vals = ();
    if ( $column_number =~ /\A (\d+) \- (\d+) \z /xms ) { 
        $val1 = $1;
        $val2 = $2;
        if ( $val1 < $val2 ) { 
            @vals = ($val1..$val2);
        }
        elsif ( $val1 > $val2 ) { 
            @vals = (reverse $val2..$val1);
        }
        elsif ( $val1 == $val2 ) { 
            warn "$column_number $val1 seems not to be an actual range.\n";
            @vals = qw( $val1 );
        }
        foreach my $valX (@vals) { 
            if ( looks_like_number($valX) ) { 
                push @human_column_numbers, $valX;
            }
            else { 
                die "$valX does not look like a number.\n";
            }
        }
    }
    elsif ( looks_like_number($column_number) ) {
        push @human_column_numbers, $column_number;
    }
    else { 
        die "Can't parse input value $column_number as a number.\n";
    }
}

# Enforce human 1-based integers:
foreach my $column_number ( @human_column_numbers ) { 
    if (    ( $column_number < 1                    ) 
         or ( $column_number != int($column_number) )
         or (! looks_like_number($column_number)    ) ) {
        die "Column number $column_number appears not to be an integer numbering from 1 to N (not from 0 or less than 0).\n";
    }
}

# Map to a Unix-style 0-based index.  Avoid weird problems with map by just copying over values.
foreach my $hum_col_val (@human_column_numbers) { 
    $hum_col_val--;
    push @perl_column_numbers, $hum_col_val;
}

while (my $input = <$INPUT_FILE>) { 
    chomp $input;
    my @unchanged_output_values = ();
    my @values_to_minmax = ();

    if ( $input !~ /\A \# /xms ) { 
        # Extract values from columns:
        my @input_values = split /\t/, $input;

        # Sanity-check the column indices:
        my $max_col_index = max(@perl_column_numbers);
        $max_col_index++;
        my $real_max_col_index = @input_values;
        if ( $real_max_col_index < $max_col_index ) { 
            die "Cannot extract $max_col_index column from data with only $real_max_col_index columns.\n";
        }

        # Convert (back) to Perl numbering:
        $max_col_index--;
        $real_max_col_index--;

        my %cols_to_minmax = ();
        foreach my $column_number (@perl_column_numbers) {
            $cols_to_minmax{$column_number} = 1;
        }

        foreach my $column_number (0..$real_max_col_index) { 
            if (! exists $cols_to_minmax{$column_number} ) { 
                push @unchanged_output_values, $input_values[$column_number];
            }
            elsif ( $cols_to_minmax{$column_number} ) { 
                push @values_to_minmax, $input_values[$column_number];
            }
            else { 
                die "Can't parse input data!\n";
            }
        }
        my $minmax_value = q{};

        if ($plain_value) { 
            if ($maximum) { 
                $minmax_value = max(@values_to_minmax);
            }
            if ($minimum) { 
                $minmax_value = min(@values_to_minmax);
            }
        }

        if ($abs_value) { 
             my @vals_sorted_by_abs = ();
             if ($maximum) {
                @vals_sorted_by_abs = sort { abs($b) <=> abs($a) } @values_to_minmax;
            }
            if ($minimum) {
                @vals_sorted_by_abs = sort { abs($a) <=> abs($b) } @values_to_minmax;
            }
            $minmax_value = $vals_sorted_by_abs[0];
        }

        if ($ratio) {
            if ($maximum) {
                $minmax_value = select_by_ratio(\@values_to_minmax, 'maximum');
            }
            if ($minimum) { 
                $minmax_value = select_by_ratio(\@values_to_minmax, 'minimum');
            }
        }
        push @unchanged_output_values, $minmax_value;
        my $output_text = join "\t", @unchanged_output_values;
        print "$output_text\n";
    }
    # Default for comment lines:
    else {
        print "$input\n";
    }
}

close $INPUT_FILE or die "Can't close filehandle to input file $infile: $!";

sub select_by_ratio { 
    my @values_to_minmax_ratio = @{ $_[0] };
    my $criterion              = $_[1];
    my @sorted_values          = ();

    if ( $criterion eq 'maximum' ) { 
        @sorted_values = sort { abs(log($b)) <=> abs(log($a)) } @values_to_minmax_ratio;
    }
    elsif ( $criterion eq 'minimum' ) { 
        @sorted_values = sort { abs(log($a)) <=> abs(log($b)) } @values_to_minmax_ratio;
    }
    else { 
        die "select_by_ratio not given correct criterion.\n";
    }
    return $sorted_values[0];
}

