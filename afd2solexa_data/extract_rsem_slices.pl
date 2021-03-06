#!/usr/bin/env perl

# extract_rsem_slices.pl -- Erich Schwarz <ems394@cornell.edu>, 2/23/2013.
# Purpose: given a single gene output file from RSEM 1.2.0, generate many single files from 1+ slices; optionally, prepend "Gene\t[TITLE]\n".

use strict;
use warnings;
use Getopt::Long;

use List::MoreUtils qw(uniq);
use Scalar::Util qw(looks_like_number);

my @initial_slice_vals = ();
my $input_file         = q{};
my $stem               = q{};
my $prefix             = q{};
my $help;

GetOptions ( 'slices=s{,}' => \@initial_slice_vals,
             'input=s'     => \$input_file,
             'output=s'    => \$stem,
             'prefix'      => \$prefix,
             'help'        => \$help,   );

if ( $help or (! @initial_slice_vals ) or (! $input_file ) ) { 
    die "Format: extract_rsem_slices.pl\n",
        "    --slices|-s  [slice value or values]\n",
        "    --input|-i   [input Acey RSEM table]\n",
        "    --output|-o  [stem name for output file]\n",
        '    --prefix|-p  [on first header line, put "Gene\t[TITLE]\n"]', "\n", 
        "    --help|-h    [print this message]\n",
        "    > [to STDOUT, Gene/value table, suitable (without header) for qual2func_table_25nov2012.pl -t]\n"; 
}
    
if (! -r $input_file) { 
    die "Can't read input file: $input_file\n";
}

my @slice_vals         = sort { $a <=> $b } @initial_slice_vals;
@slice_vals            = uniq @slice_vals;

foreach my $slice_val (@slice_vals) { 
    # Reality-check the numbers.
    if (! looks_like_number($slice_val) ) { 
        die "slice_val $slice_val does not look like number\n";
    }
    if ( $slice_val != int($slice_val) ) { 
        die "slice_val $slice_val does not look like integer\n";
    }
    if ( $slice_val <= 1 ) { 
        die "slice_val $slice_val is less than 2\n";
    }

    my @slice_inputs = `cut -f 1,$slice_val $input_file`;
    chomp @slice_inputs;
    my $header = shift @slice_inputs;
    if ( $header =~ /\A (?:Gene|gene_id) \t (\S+) \z/xms ) { 
        my $title = $1;
        if ($prefix) {
            $prefix = "Gene\t$stem.$title";
        }

        # Alter special characters so that they won't cripple a filename:
        $title =~ s/\//.vs./g;
        $title =~ s/\W/./g;
        # Get rid of redundant '.'s:
        $title =~ s/[.]+/./g;

        my $output_file = $stem . q{.} . $title . '.txt';
        # Again, get rid of redundant '.'s:
        $output_file =~ s/[.]+/./g;
        # Enforce no overwriting:
        $output_file = safename($output_file);

        open my $OUT, '>', $output_file or die "Can't open output file $output_file\n";
        print $OUT "$prefix\n" if $prefix;
        foreach my $slice_input_line (@slice_inputs) {
            print $OUT "$slice_input_line\n";
        }
        close $OUT or die "Can't close filehandle to output file $output_file\n";
    } 
    else { 
        die "Can't parse header of slice inputs: $header\n";
    }
}


sub safename {
    my $filename = $_[0];
    my $orig_filename = $filename;
    if (-e $orig_filename) {
        my $suffix1 = 1;
        $filename = $filename . ".$suffix1";
        while (-e $filename) {
            $suffix1++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix1";
        }
    }
    return $filename;
}

