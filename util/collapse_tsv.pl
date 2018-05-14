#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Getopt::Long;
use List::MoreUtils qw(uniq);

my $infile = q{};
my $help;

my @key_items  = ();
my @data_types = ();

my $header_count = 0;
my $data_count   = 0;

my $header = q{};

my $data_ref;

GetOptions ( 'infiles=s' => \$infile,
             'help'         => \$help,   );

if ( $help or (! $infile) ) { 
    die "Format: count_fasta_residues.pl\n",
        "    --infile|-i   <input stream/file>\n",
        "    --help|-h     [print this message]\n",
        ;
}

my $INPUT_FILE;
if ($infile eq '-') {
    # Special case: get the stdin handle
    $INPUT_FILE = *STDIN{IO};
}
else {
    # Standard case: open the file
    open $INPUT_FILE, '<', $infile or die "Can't open input file $infile. $!\n";
}

while (my $input = <$INPUT_FILE>) { 
    chomp $input;
    if (! $header ) {
        $header = $input;
        @data_types = split /\t/, $input;

        my %seen = ();
        foreach my $column (@data_types) {
            if ( exists $seen{$column} ) { 
                die "Redundant column name (\"$column\") in header line: $input\n";
            }
            $seen{$column} = $1;
            if ( $column !~ / \S /xms ) { 
                die "Column name (\"$column\") needs at least one non-space character in header line: $input\n";
            }
        }

        $header_count = @data_types;
    }
    else {
        my @data_entries = split /\t/, $input;
        $data_count = @data_entries;
        if ( $data_count != $header_count ) {
            die "Data line has $data_count entries, but header has $header_count entries -- offending data is: $input\n";
        }

        $data_count--;

        my $key_item = $data_entries[0];
        push @key_items, $key_item;

        # start with $i *1*, not *0*:
        foreach my $i (1..$data_count) { 
            my $data_type  = $data_types[$i];
            my $data_entry = $data_entries[$i];
            $data_ref->{'key_item'}->{$key_item}->{$data_type}->{$data_entry} = 1;
        }
    }
}
close $INPUT_FILE;

@key_items = uniq(@key_items);

$header_count--;

foreach my $key_item (@key_items) {
    print "$header\n" if $header;
    $header = q{};

    my @data_entries = ();
    push @data_entries, $key_item;

    foreach my $i (1..$header_count) {
        my $data_type  = $data_types[$i];
        my @data_type_entries = sort keys %{ $data_ref->{'key_item'}->{$key_item}->{$data_type} };
        my $data_entry = join "; ", @data_type_entries;
        push @data_entries, $data_entry;
    }

    my $output = join "\t", @data_entries;
    print "$output\n";
}

