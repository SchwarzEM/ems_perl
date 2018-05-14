#!/usr/bin/env perl

# shift_preExcel_data_column.pl -- Erich Schwarz <emsch@its.caltech.edu>, 7/24/2011.
# Purpose: given a list of genes to shift, take a pre-Excel TSV file and either add an empty extra column (with "\t") or -- for gene in list -- move rightmost data column to newly created column on right edge of table in rightmost column ("[added \t] [data line]").

use strict;
use warnings;

my $gene     = q{};
my $mid_text = q{};
my $data_col = q{};
my $end_text = q{};

my $genelist  = $ARGV[0];
my $data_file = $ARGV[1];

my %genes2shift = ();

open my $GENES, '<', $genelist or die "Can't open gene list $genelist: $!";
while (my $input = <$GENES>) { 
    chomp $input;
    if ( $input =~ /\A (\S+) \z /xms ) { 
        $gene = $1;
        $genes2shift{$gene} = 1;
    }
}
close $GENES or die "Can't close filehandle to gene list $genelist: $!";

open my $DATA, '<', $data_file or die "Can't open data file $data_file: $!";
while (my $input = <$DATA>) { 
    chomp $input;
    my $output = q{};
    if ( $input =~ / \A (\S+) \t (.*?) (\S+) (\s*) \z /xms ) { 
        $gene     = $1;
        $mid_text = $2;
        $data_col = $3;
        $end_text = $4;
        if (exists $genes2shift{$gene}) { 
            $output = $gene . "\t" . $mid_text . $end_text . "\t" . $data_col;
        }
        else { 
            $output = $input . "\t";
        }
        print "$output\n";
    }
    else { 
        warn "Can't parse input line: $input\n";
        print "$input\n";
    }
}

