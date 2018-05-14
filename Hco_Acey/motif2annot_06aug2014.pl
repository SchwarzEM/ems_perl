#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $motif_list = $ARGV[0];
my $pfam_table = $ARGV[1];

my %seen = ();

open my $MOTIF, '<', $motif_list;
while (my $input = <$MOTIF>) { 
    chomp $input;
    if ( $input =~ /\A \S+ \z/xms ) {
        $seen{$input} = 1;
    }
    else {
        die "From motif list $motif_list, cannot parse: $input\n";
    }
}
close $MOTIF;

open my $PFAM, '<', $pfam_table;
while (my $input = <$PFAM>) {
    chomp $input;
    if ( $input =~ /\A \S+ \s+ \S+ \s+ \S+ \s+ (\S+) /xms ) {
        my $motif   = $1;
        if ( $seen{$motif} ) { 
            print "$input\n";
        }
    }
    elsif ( $input !~ /\A [#] /xms ) { 
        die "From PFAM hit table $pfam_table, cannot parse: $input\n";
    }
}
close $PFAM;

