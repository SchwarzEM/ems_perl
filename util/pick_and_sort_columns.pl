#!/usr/bin/env perl

# pick_and_sort_columns.pl -- Erich Schwarz <emsch@its.caltech.edu>, 6/5/2011.
# Purpose: from a tab-defined input file or stream, pick and print 1+ named tab-defined columns; pass on comment ('# ...') lines untouched; significant fix on 6/5/2011 to deal with Perl's inability to Do The Right Thing on splitting a tab-delimited line with empty columns.

use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);
use List::Util qw(max);

my $infile               = q{};
my @input_column_numbers = ();
my @human_column_numbers = ();
my @perl_column_numbers  = ();
my $help;

GetOptions ( 'column_numbers=s{,}' => \@input_column_numbers,
             'infile=s'            => \$infile,
             'help'                => \$help, );

if ( $help or (! $infile ) or (! @input_column_numbers ) ) { 
    die "Format: pick_and_sort_columns.pl --infile|-i <input stream/files>  --column_numbers|-c [columns of input to print, IN ORDER LISTED] --help|-h\n";
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

# Interpret the *strings* given to --column_numbers, so that ranges like '1-4' or '10-7' can be interpreted.
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
            warn "$column_number seems not to be an actual range.\n";
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
    my @input_values  = ();
    my @output_values = ();

    if ( $input !~ /\A \# /xms ) { 
        # Extract values from column in way that doesn't treat empty fields as uninitialized, since Perl doesn't manage it with simple "split "\t", $input"!
        my $input_copy     = $input;
        my $new_value      = q{};
        my $remaining_line = q{};


        # Enforce deterministic parsing here, and two other places shortly below.
        # Note that this parse doesn't *require* a lot, but it does demand that the parser accept this regex as OK:        
        if ( $input_copy !~ /\A [^\t]* .*? \z/xms ) {
            die "Can't parse: $input_copy\n";
        }
        if ( $input_copy =~ /\A ([^\t]*) (.*?) \z/xms ) { 
            $new_value      = $1;
            $remaining_line = $2;
            push @input_values, $new_value;
        }

        # Note that this parse requires that the input line should have had at least *two* columns separated by at least one tab.
        #     That is reasonable since we are not likely to use this script on single-column data, but it is a slight loss of elegance and generality.
        #     (If only, if only Perl coped with splitting the line on "\t"!)
        if ( $remaining_line !~ /\A \t [^\t]* /xms ) {
            die "In input line (\"$input_copy\"), can't parse remainder: \"$remaining_line\"\n";
        }

        # For two or more leading tabs, step through each non-tab field following a tab (even if the field's empty space) and explicitly populate @input_values.
        #     $new_value has been initialized earlier, so it should remain an existing variable (though empty sometimes).
        #     Note that this should get skipped over completely if there is just one tab and one non-tab field left.
        while ($remaining_line =~ /\A \t ([^\t]*) (\t .*) /xms ) { 
            $new_value      = $1;
            $remaining_line = $2;
            push @input_values, $new_value;
        }

        # Eventually this should pare down $remaining_line to tab-less singularity.  Enforce that, then import it.
        if ( $remaining_line !~ /\A \t [^\t]* \z /xms ) {
            die "In input line (\"$input_copy\"), can't parse remainder: \"$remaining_line\"\n";
        }
        if ( $remaining_line =~ /\A \t ([^\t]*) \z /xms ) {
            $new_value      = $1;
            push @input_values, $new_value;
        }

        # Sanity-check the column indices:
        my $max_col_index = max(@perl_column_numbers);
        $max_col_index++;

        # Old method for counting columns in data, which turns out to fail if data columns are empty:
        my $orig_real_max_col_index = @input_values;

        # New method for counting columns in data, which should be bulletproof to empty data columns:
        my $real_max_col_index = ($input =~ tr/\t/\t/);
        # But, note that N columns will be split by N-1 tabs!  So, correct accordingly:
        $real_max_col_index++;

        # Error-checking to make sure I am properly populating @input_values:
        if ( $real_max_col_index != $orig_real_max_col_index ) {
            warn "Getting discrepancy between seeing $real_max_col_index columns with new method and $orig_real_max_col_index with old one.\n";
        }

        # Fail loudly if the user asks for nonexistent columns:
        if ( $real_max_col_index < $max_col_index ) { 
            die "Cannot extract $max_col_index column from data with only $real_max_col_index columns.\n";
        }

        # At last!  The point of the whole exercise!
        foreach my $column_number (@perl_column_numbers) {
            push @output_values, $input_values[$column_number];
        }

        my $output_text = join "\t", @output_values;
        print "$output_text\n";
    }
    # Default for comment lines:
    else {
        print "$input\n";
    }
}

close $INPUT_FILE or die "Can't close filehandle to input file $infile: $!";

