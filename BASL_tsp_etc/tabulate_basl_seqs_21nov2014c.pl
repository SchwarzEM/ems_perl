#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data     = $ARGV[0];
my $nlses    = $ARGV[1];
my %seen_nls = ();

open my $NLSES, '<', $nlses;
while (my $input = <$NLSES>) { 
    chomp $input;
    if ( $input =~ /\A (\S+) \s+ posterior \s+ /xms ) {
        my $seq = $1;
        $seen_nls{$seq} = 1;
    }
}
close $NLSES;

open my $DATA, '<', $data;
while (my $input = <$DATA>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \s .+ \z/xms ) { 
        my $seq = $1;
        my $added_annot = "\t";
        if ( exists $seen_nls{$seq} ) {
            $added_annot = "\tPredicted_NLS";
        }
        print "$input$added_annot\n";
    }
    else {
        die "From data table $data, cannot parse input: $input\n";
    }
}
close $DATA;

