#!/usr/bin/env perl

# emp_pzero_cols.pl -- Erich Schwarz <emsch@its.caltech.edu>, 8/12/2012.
# Purpose; given a tab-delimited file with numbers in columns, map all zero values in a column to the smallest observed non-zero value in that column.  Be careful that data from different organisms, tissues, etc., are not mixed into one column before doing this!

use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my $data_ref;

my $infile = q{};
my $help;

GetOptions ( 'infile=s' => \$infile,
             'help'     => \$help,   );

if ( $help or (! $infile ) ) {
die "Format: emp_pzero_cols.pl\n",
    "    --infile|-i   <single input file>\n",
    "    --help|-h     [print this message]\n",
    ;
}

open my $INPUT_FILE, '<', $infile or die "Can't open input file $infile. $!\n";

while (my $input = <$INPUT_FILE>) { 
    chomp $input;
    my @input_vals = split /\t/, $input;
    my $input_index = @input_vals;
    $input_index--;
    foreach my $i (0..$input_index) {
        if ( ( looks_like_number($input_vals[$i]) ) and ( $input_vals[$i] > 0 ) ) { 
            if ( ( exists $data_ref->{'column'}->{$i}->{'pzero'} ) and ( $data_ref->{'column'}->{$i}->{'pzero'} > $input_vals[$i] ) ) { 
                $data_ref->{'column'}->{$i}->{'pzero'} = $input_vals[$i];
            }
            if (! exists $data_ref->{'column'}->{$i}->{'pzero'} ) { 
                $data_ref->{'column'}->{$i}->{'pzero'} = $input_vals[$i];
            }
       }
   }
}

close $INPUT_FILE or die "Can't close filehandle to input file $infile. $!\n";

open $INPUT_FILE, '<', $infile or die "Can't open input file $infile. $!\n";

while (my $input = <$INPUT_FILE>) {
    chomp $input;
    my @input_vals = split "\t", $input;
    my $input_index = @input_vals;
    $input_index--;
    foreach my $i (0..$input_index) {
        if ( ( looks_like_number($input_vals[$i]) ) and ( $input_vals[$i] == 0 ) ) {
            $input_vals[$i] = $data_ref->{'column'}->{$i}->{'pzero'};
        }
    }
    my $output = join "\t", @input_vals;
    print "$output\n";
}   

close $INPUT_FILE or die "Can't close filehandle to input file $infile. $!\n";

