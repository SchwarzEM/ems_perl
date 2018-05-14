#!/usr/bin/env perl

# aceify_vert_mot.pl -- Erich Schwarz <emsch@its.caltech.edu>, 4/1/2010.
# Purpose: take very basic vertical PFM/PWM text and print out the core of a .ace-able motif; dynamically assign residues to data columns, so non-standard input orders (e.g., "A T C G") work OK.

use strict;
use warnings;
use Getopt::Long;

my $site_vals_ref;
my $help;

my $ok_res = 'aAcCgGtTuU';

my $motif_no         = 0;
my $read_nos         = 0;
my $filler_columns   = 0;
my $residues_subline = q{};
my @column_residues  = ();

GetOptions ( 'filler_columns=i' => \$filler_columns,
             'help'             => \$help,           );

my $fill_plus_cols = $filler_columns + 3;

if ($help) { 
    die "Format: aceify_vert_mot.pl",
        " --filler_columns|-f [no. of right-hand columns to skip]",
        " --help|-h <input stream/file>",
        "\n",
        ;
}

while (my $input = <>) { 
    chomp $input;
    if ( $read_nos  
         and ( $input =~ / \A 
                           (?:\S+\s+){$filler_columns} 
                           (\S+) 
                           \s+ 
                           (\S+) 
                           \s+ 
                           (\S+)
                           \s+ 
                           (\S+) 
                           \s* /xms ) ) {
        my $val_col_0 = $1;
        my $val_col_1 = $2; 
        my $val_col_2 = $3;
        my $val_col_3 = $4;
        push @{ $site_vals_ref->{$motif_no}->{$column_residues[0]} }, $val_col_0;
        push @{ $site_vals_ref->{$motif_no}->{$column_residues[1]} }, $val_col_1;
        push @{ $site_vals_ref->{$motif_no}->{$column_residues[2]} }, $val_col_2;
        push @{ $site_vals_ref->{$motif_no}->{$column_residues[3]} }, $val_col_3;
    }
    if ( $read_nos 
         and ( $input !~ / \A
                           (?:\S+\s+){$fill_plus_cols}
                           (\S+) \s* /xms ) ) {
        $read_nos = 0;
    }
    if ( $input =~ /\A (?:\S+\s+){$filler_columns} 
                       ( [$ok_res] 
                         \s+ 
                         [$ok_res] 
                         \s+ 
                         [$ok_res] 
                         \s+ 
                         [$ok_res] ) 
                   /xms ) { 
        $residues_subline = $1;
        @column_residues = split /\s+/, $residues_subline;
        $residues_subline = join q{}, (sort (split /\s+/, $residues_subline) );
        $residues_subline = uc $residues_subline;
        if ( ( $residues_subline eq 'ACGT' ) or ( $residues_subline eq 'ACGU' ) ) { 
            $read_nos = 1;
            $motif_no++;
        }
    }
}

foreach my $motif ( sort { $a <=> $b } keys %{ $site_vals_ref } ) { 
    print "\n", q{Position_Matrix : "WBPmat00000xxx_}, $motif, q{"}, "\n";
    foreach my $res ( qw( A C G T ) ) {
        if ( exists $site_vals_ref->{$motif}->{$res} ) { 
            my $value_line = join q{  }, @{ $site_vals_ref->{$motif}->{$res} };
            print 'Site_values    ', "$res", q{    }, "$value_line\n";  
        }
    }
    print "\n";
}

