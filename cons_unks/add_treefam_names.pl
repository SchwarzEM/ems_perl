#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $treefam_names = q{};
my $treefam_table = q{};

$treefam_names = $ARGV[0] if $ARGV[0];
$treefam_table = $ARGV[1] if $ARGV[1];

my $data_ref;

if ( (! -e $treefam_names ) or (! -e $treefam_table ) ) {
    die "Format: add_treefam_names.pl [table with TreeFam IDs with human-readable names] [table with TreeFam IDs lacking names] > [2cd table, with names added]\n";
}

open my $NAMES, '<', $treefam_names;
while (my $input = <$NAMES>) {
    chomp $input;

    # Sample input lines:
    # TF_name	fam_symbol	fam_description [etc.]
    # TF333215	NaN	NaN [etc.]
    # TF105964	NP_055412.2	estrogen receptor binding protein [etc.]

    if ( $input =~ /\A (TF\d+) \t [^\t]* \t ([^\t]+) \t /xms ) { 
        my $treefam_id    = $1;
        my $treefam_annot = $2;
        if ( $treefam_annot ne 'NaN' ) {
            $treefam_annot =~ s/[ ]+/ /g;
            $treefam_annot =~ s/\A\s+//;
            $treefam_annot =~ s/\s+\z//;
            $treefam_annot = q{"} . $treefam_annot . q{"};
            $data_ref->{'treefam_id'}->{$treefam_id}->{'treefam_annot'} = $treefam_annot;
        }
    }
}
close $NAMES;

open my $TABLE, '<', $treefam_table;
while (my $input = <$TABLE>) {
    chomp $input;
    if ( $input !~ / TF\d+ /xms ) {
        print "$input\n";
    }
    else {
        my $leading_text  = q{};
        my $treefam_id    = q{};
        my $trailing_text = q{};
        my $output        = q{};

        # Do dynamic and conditional substitutions, down the line of the text:
        while ( $input =~ /\A (.*?) (TF\d+) (.*) \z/xms ) { 
            $leading_text  = $1;
            $treefam_id    = $2;
            $trailing_text = $3;

            # Only add a text annotation *if* it really exists.
            if ( exists $data_ref->{'treefam_id'}->{$treefam_id}->{'treefam_annot'} ) { 
                $treefam_id = $treefam_id . q{|} . $data_ref->{'treefam_id'}->{$treefam_id}->{'treefam_annot'};
            }
            $output .= $leading_text;
            $output .= $treefam_id;
            $input  = $trailing_text;
        }

        # Once the substitutions are over:
        $output .= $input;
        print "$output\n";
    }
}
close $TABLE;


