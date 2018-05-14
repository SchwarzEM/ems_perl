#!/usr/bin/env perl

# shift_preExcel_data_column.pl -- Erich Schwarz <emsch@its.caltech.edu>, 7/10/2012.
# Purpose: given a list of genes to shift, take a pre-Excel TSV file and either add an empty extra column (with "\t") or -- for gene in list -- move rightmost data column to newly created column on right edge of table in rightmost column ("[added \t] [data line]").

use strict;
use warnings;
use Getopt::Long;

my $data_file       = q{};
my @genelists       = ();
my $only_unshifted;
my $help;

my $gene            = q{};
my $x_axis          = q{};
my $y_axis          = q{};
my $end_tabs        = q{};

my @orig_data_lines = ();
my @rev_data_lines  = ();

my %genes2shift     = ();

GetOptions ( 'data=s'     => \$data_file,
             'genes=s{,}' => \@genelists,
             'only'       => \$only_unshifted,
             'help'       => \$help, );

if ( $help or (! @genelists ) or (! $data_file ) ) { 
    die "Format: shift_preExcel_data_column.pl\n",
        "            --data|-d <data_file to recolumn>\n",
        "            --genes|-g <gene_lists -- can have 1, 2, or more; each one gets indented in turn>\n",
        "            --only|-o [only touch previously unshifted columns]\n",
        "            --help|-h \n",
        ;
}

# Before bothering with anything else, do two things.
#
#     Ensure that the original data file has a single width of tabs:
#     Import the original data file to a stored array, which can then be indented 2+ times.

open my $DATA, '<', $data_file or die "Can't open data file $data_file: $!";
my $tab_count      = 0;
my $prev_tab_count = 0 ;

while (my $input = <$DATA>) {
    chomp $input;
    if ( $tab_count > 0 ) {
        $prev_tab_count = $tab_count;
    }
    $tab_count = ( $input =~ tr/\t/\t/ );

    # Enforce equal tabcount throughout data file:
    if ( ( $prev_tab_count > 0 ) and ( $prev_tab_count != $tab_count ) ) { 
        die "Data file $data_file has inconsistent numbers of tabs-per-line ($prev_tab_count versus $tab_count) at input line: $input\n";
    }
    # Archive original data file's lines:
    push @orig_data_lines, $input;
}
close $DATA or die "Can't close filehandle to data file $data_file: $!";

# For each list of genes, carry out an indentation in turn.
foreach my $genelist (@genelists) { 
    if (@rev_data_lines) { 
        @orig_data_lines = @rev_data_lines;
    }
    @rev_data_lines  = ();

    # For *each* genelist, this hash must be rebuilt from zero.
    %genes2shift = ();

    open my $GENES, '<', $genelist or die "Can't open gene list $genelist: $!";

    while (my $input = <$GENES>) { 
        chomp $input;
        if ( $input =~ /\A (\S+) \z /xms ) { 
            $gene = $1;
            $genes2shift{$gene} = 1;
        }
        else { 
            die "From gene list $genelist, can't parse input line: $input\n";
        }
    }
    close $GENES or die "Can't close filehandle to gene list $genelist: $!";

    foreach my $input (@orig_data_lines) { 
        chomp $input;
        my $output = q{};
        if (! $only_unshifted ) { 
            if ( $input =~ / \A (\S+) \t (\S+ \t+) (\S+) (\t*) \z /xms ) { 
                $gene     = $1;
                $x_axis   = $2;
                $y_axis   = $3;
                $end_tabs = $4;
                if (exists $genes2shift{$gene}) { 
                    $output = $gene . "\t" . $x_axis . $end_tabs . "\t" . $y_axis;
                }
                else { 
                    $output = $gene . "\t" . $x_axis . $y_axis . $end_tabs . "\t";
                }
                push @rev_data_lines, $output;
            }
            else {
                die "From data file $data_file, can't parse input line: $input\n";
            }
        }
        elsif ( $only_unshifted ) {
            if ( $input =~ / \A (\S+) \t (\S+ \t) (\S+) (\t*) \z /xms ) {
                $gene     = $1;   
                $x_axis   = $2;
                $y_axis   = $3;
                $end_tabs = $4;
                if (exists $genes2shift{$gene}) {
                    $output = $gene . "\t" . $x_axis . $end_tabs . "\t" . $y_axis;
                }
                else {
                    $output = $gene . "\t" . $x_axis . $y_axis . $end_tabs . "\t";
                }
                push @rev_data_lines, $output;
            }
            elsif ( $input =~ / \A \S+ \t \S+ \t{2,} \S+ \t* \z /xms ) {
                $output = $input . "\t";
                push @rev_data_lines, $output;
            }
            else {
                die "From data file $data_file, can't parse input line: $input\n";
            }
        }
        else {
            die "Unexplained failure to define choice of shifting mode.\n";
        }
    }
}

foreach my $output (@rev_data_lines) { 
    print "$output\n";
}

