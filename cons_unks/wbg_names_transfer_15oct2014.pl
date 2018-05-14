#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

# wbg_names_transfer_15oct2014.pl -- Erich Schwarz <ems394@cornell.edu>, 10/15/2014.
# Purpose: given new WBGene names in one file and a table that has some kind of WBGene names in another, use the former to update the latter's names.

my $new_names_file = q{};
my $current_table  = q{};

$new_names_file = $ARGV[0] if $ARGV[0];
$current_table  = $ARGV[1] if $ARGV[1];

if ( (! -r $new_names_file) or (! -r $current_table ) ) {
    die "Format: wbg_names_transfer_15oct2014.pl [new names file] [current table] > [name-updated table]\n";
}

my $data_ref;

open my $NEW_NAMES, '<', $new_names_file;
while (my $input = <$NEW_NAMES>) {
    chomp $input;

    # Sample input lines:
    # WBGene00000001	aap-1	Y110A7A.10
    # WBGene00000308	C04A11.t1	C04A11.t1
    if ( $input =~ /(WBGene\d+)\t(\S+)\t(\S+)/xms ) { 
        my $wbgene_id = $1;
        my $cgc_id    = $2;
        my $cds_id    = $3;
        my $new_name  = q{};
        if ( ( exists $data_ref->{'wbgene_id'}->{$wbgene_id}->{'new_name'} ) and ( $data_ref->{'wbgene_id'}->{$wbgene_id}->{'new_name'} ne $new_name ) ) {
            die "In new names table $new_names_file, inconsistent assignment of WB Gene ID $wbgene_id",
                " to both $data_ref->{'wbgene_id'}->{$wbgene_id}->{'new_name'} and $new_name\n",
            ;
        }
        $new_name = $wbgene_id . q{|} . $cds_id;
        if ( $cgc_id ne $cds_id ) {
            $new_name = $new_name . q{|} . $cgc_id;
        }
        $data_ref->{'wbgene_id'}->{$wbgene_id}->{'new_name'} = $new_name;
    }
    else { 
        die "From new names table $new_names_file, cannot parse: $input\n";
    }
}
close $NEW_NAMES;

open my $CURRENT_TABLE, '<', $current_table;
while (my $input = <$CURRENT_TABLE>) {
    chomp $input;

    # To avoid making the replacement repeatedly trip over itself by replacing a shorter name with a longer name (!),
    #    compel it to work down the line of text; do not allow it to see the same text twice.
    # For each line of the table, allow global replacements to happen -- but only one time per new name, to avoid an infinite loop.

    my $output = q{};

    while ( $input =~ /\A (.*?) ((WBGene\d+)\S+) \b (.*) \z/xmsg ) { 
        my $leading_text   = $1;
        my $old_name_text  = $2;   # Distinguish this from the old *name*; we want to alter the text, but keep unchanged memory of what the old name was.
        my $wbgene_id      = $3;
        my $following_text = $4;
        my $new_name       = q{};

        my $old_name = $old_name_text;

        if ( exists $data_ref->{'wbgene_id'}->{$wbgene_id}->{'new_name'} ) {
            $new_name = $data_ref->{'wbgene_id'}->{$wbgene_id}->{'new_name'};
            if ( $new_name ne $old_name ) {
                # Block '|' characters from being misinterpreted as pattern elements!
                $old_name =~ s/[|]/[\][|]/g;
                $old_name =~ s/\[\]//g;   # For some idiotic reason, Perl is putting '[]' into my patternized text.

                $old_name_text =~ s/$old_name/$new_name/;
            }
        }
        else {
            warn "Was unable to identify possible new name for $old_name\n";
        }

        # Append the leading text and the now-substituted old name text to what will become our output.
        $output = $output . $leading_text . $old_name_text ;

        # Move things forward by redefining $input as our *remaining* text.
        $input = $following_text;
    }

    # Add back any unprocessed $input to $output (which can be *all* of original $input, if zero substitutions happened).
    $output = $output . $input;
    print "$output\n";
}

