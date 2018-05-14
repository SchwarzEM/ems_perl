#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data      = $ARGV[0];
my $basls     = $ARGV[1];
my %seen_basl = ();

open my $BASLS, '<', $basls;
while (my $input = <$BASLS>) { 
    chomp $input;
    $input =~ s/\s//g;
    $seen_basl{$input} = 1;
}
close $BASLS;

open my $DATA, '<', $data;
while (my $input = <$DATA>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t .+ \z/xms ) { 
        my $seq = $1;
        my $added_annot = "\t";
        if ( exists $seen_basl{$seq} ) {
            $added_annot = "\tBASL_subfamily";
        }
        print "$input$added_annot\n";
    }
    else {
        die "From data table $data, cannot parse input: $input\n";
    }
}
close $DATA;

