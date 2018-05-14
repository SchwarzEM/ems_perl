#!/usr/bin/env perl

# append_clean_data_to_sloppy_table.pl -- Erich Schwarz <ems394@cornell.edu>, 7/31/2015.
# Purpose: given a clean Gene annotation table, append its data to a sloppy (redundant-lined) Gene annotation table.

use strict;
use warnings;
use autodie;

my $sloppy_table = q{};
my $clean_table  = q{};

$sloppy_table = $ARGV[0] if $ARGV[0];
$clean_table  = $ARGV[1] if $ARGV[1];

if ( (! -r $sloppy_table) or (! -r $clean_table) ) {
    die "Format: append_clean_data_to_sloppy_table.pl [sloppy primary table] [clean appendable table] > [STDOUT or a file]\n";
}

my %header_ok  = ();
my %gene2annot = ();

$header_ok{$sloppy_table} = 'maybe';
$header_ok{$clean_table}  = 'maybe';

my $appended_tabs = 0;

open my $CLEAN, '<', $clean_table;
while (my $input = <$CLEAN>) {
    chomp $input;
    if ( $input =~ /\A (\S+) ( \t [^\t]* (?: \t [^\t]*)* ) \z/xms ) { 
        # The very first of these needs to be 'Gene', so that we have a standard header line.
        my $gene_id    = $1;

        # This need not have any actual text, but there *must* be at least one \t after $gene_id.
        my $gene_annot = $2;

        # Enforce there being a top Gene line.
        if ( $header_ok{$clean_table} eq 'maybe' ) {
            if ( $gene_id ne 'Gene' ) {
                die "Very first line in $clean_table needs to start with \"Gene\", not like this: \"$input\"\n";
            }
            else {
                $header_ok{$clean_table} = 'yes';
            }
        }

        # We tolerate sloppiness in the primary file, but *not* the appended file.
        if ( exists $gene2annot{$gene_id} ) {
            die "In allegedly clean added annotation file $clean_table, redundant annotation for $gene_id: \"$input\"\n";
        }

        # Each appended annotation has to have exactly the same number of tabs; and we need to know *what* that number is,
        #     so that (if necessary) we can create an empty annotation line consisting of nothing but tabs!
        my $input_tab_count = ( $gene_annot =~ tr/\t/\t/ );
        if ( ( $appended_tabs > 0 ) and ( $appended_tabs != $input_tab_count ) ) {
            die "Contradictory non-zero tab counts in appended annotation: $appended_tabs versus $input_tab_count\n";
        }
        else {
            # OK to repeat this with every single line, as long as the numbers stay unchanged!
            $appended_tabs = $input_tab_count;
        }

        # Having passed all the tests, record the annotation for later appendage.
        $gene2annot{$gene_id} = $gene_annot;
    }
    # Enforce successful parsing of input lines:
    else {
        die "In allegedly clean added annotation file $clean_table, cannot parse input line: \"$input\"\n";
    }
}
close $CLEAN;

# Generate a blank annotation line that we can add to genes that don't have any annotation in our appended file.
my $default_annotation = "\t" x $appended_tabs;

open my $SLOPPY, '<', $sloppy_table;
while (my $input = <$SLOPPY>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t [^\t]* (?: \t [^\t]* )* \z/xms ) {
        # The very first of these needs to be 'Gene', so that we have a standard header line.
        my $gene_id    = $1;
        
        # Enforce there being a top Gene line.
        if ( $header_ok{$sloppy_table} eq 'maybe' ) {
            if ( $gene_id ne 'Gene' ) {
                die "Very first line in $sloppy_table needs to start with \"Gene\", not like this: \"$input\"\n";
            }
            else {
                $header_ok{$sloppy_table} = 'yes';
            }
        }

        # By default, we add nothing; but if there *is* available annotation, then we add that instead.
        my $appended_annotation = $default_annotation;
        if ( exists $gene2annot{$gene_id} ) {
            $appended_annotation = $gene2annot{$gene_id}
        }
        
        # Finally!
        print "$input\t$appended_annotation\n";
    }   
    # Enforce successful parsing of input lines:
    else {
        die "In primary (sloppy) annotation file $sloppy_table, cannot parse input line: \"$input\"\n";
    }   
}       
close $SLOPPY;

