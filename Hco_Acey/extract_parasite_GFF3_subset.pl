#!/usr/bin/env perl

# extract_parasite_GFF3_subset.pl -- Erich Schwarz <ems394@cornell.edu>, 11/14/2022.
# Purpose: w/ listed seqs. or 1-name argument, extract subset of ParaSite GFF3; warn of misses.

use strict;
use warnings;
use autodie;

unless ($#ARGV == 1) { 
    die "Format: extract_parasite_GFF3_subset.pl",
        "  [gene list file OR gene name]",
        "  [large ParaSite GFF3 to extract]\n",
        ;
}

my $input_list     = $ARGV[0];
my $input_gff3     = $ARGV[1];
my %input_names    = ();
my $reading_subset = "no";
my $warnings       = $input_list . ".warnings";

if (-e $input_list) { 
    open (my $INPUT_LIST, "$input_list");

    while (my $inline = <$INPUT_LIST>) { 
        chomp $inline;
        if ( $inline =~ /\A (\S+) \z/xms) { 
            my $gene = $1;
            $input_names{$gene} = 'wanted';
        }
        else {
            die "In input list $input_list, cannot parse input line: $inline\n";
        }
    }
    close $INPUT_LIST;
}

# If no list file, treat as single-name argument.
if (! -e $input_list) { 
    $input_names{$input_list} = 'wanted';
}

open my $INPUT_GFF3, '<', $input_gff3;

while (my $gff3_line = <$INPUT_GFF3>) {
    if ( ( $gff3_line =~ /\A [#][ ]Gene[ ]gene[:](\S+) /xms ) and ( $input_names{$1} ) ) { 
        my $gene = $1;
        $input_names{$gene} = 'found';        
        print "$gff3_line";
        $reading_subset = 'yes';
    }
    elsif ( ( $gff3_line =~ /\A [#][ ]Gene[ ]gene[:](\S+) /xms ) and (! $input_names{$1} ) ) {
        $reading_subset = 'no';
    }
    elsif ( $reading_subset eq 'yes' ) {
        print "$gff3_line";
    }
}
close $INPUT_GFF3;

# Scan hash for any non-zero values; warn that these were missed.

open my $WARNINGS, '>', "$warnings";
foreach my $key (sort keys %input_names) {
    if ( $input_names{$key} eq 'wanted' ) { 
        print {$WARNINGS} "$key not found in $input_gff3\n";
    }
}
close $WARNINGS;
