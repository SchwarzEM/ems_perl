#!/usr/bin/env perl

# shift_preExcel_data_column_24jul2011.pl -- Erich Schwarz <emsch@its.caltech.edu>, 7/24/2011.
# Purpose: [LEGACY VERSION] given a list of genes to shift, take a pre-Excel TSV file and either add an empty extra column (with "\t") or -- for gene in list -- move rightmost data column to newly created column on right edge of table in rightmost column ("[added \t] [data line]").

use strict;
use warnings;
use Getopt::Long;

my $gene        = q{};
my $mid_text    = q{};
my $data_col    = q{};
my $end_text    = q{};
my $genelist    = q{};
my $data_file   = q{};
my $x_axis      = q{};
my $y_axis      = q{};
my %genes2shift = ();

my $only_unshifted;

my $help;

GetOptions ( 'genes=s' => \$genelist,
             'data=s'  => \$data_file,
             'only'    => \$only_unshifted,
             'help'    => \$help, );

if ( $help or (! $genelist ) or (! $data_file ) ) { 
    die "Format: shift_preExcel_data_column.pl\n",
        "            --genes|-g <gene_list>\n",
        "            --data|-d <data_file to recolumn>\n",
        "            --only|-o [only touch previously unshifted columns]\n",
        "            --help|-h \n",
        ;
}

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
    if (! $only_unshifted ) { 
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
    elsif ( $only_unshifted ) {
        if ( $input =~ / \A (\S+) \t (S+) \t (\S+) (.*) \z /xms ) {
            $gene     = $1;   
            $x_axis   = $2;
            $y_axis   = $3;
            $end_text = $4;
            if (exists $genes2shift{$gene}) {
                $output = $gene . "\t" . $x_axis . "\t" . $end_text . "\t" . $y_axis;
            }
            else {
                $output = $input . "\t";
            }
            print "$output\n";
        }
        elsif ( $input =~ / \A \S+ \t .*? \S+ \s* \z /xms ) {
            $output = $input . "\t";
            print "$output\n";
        }
        else {
            warn "Can't parse input line: $input\n";
            print "$input\n";
        }
    }
    else {
        die "Unexplained failure to define choice of shifting mode.\n";
    }
}

