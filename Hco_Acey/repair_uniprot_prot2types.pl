#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $full_names = $ARGV[0];
my $uprot_dat  = $ARGV[1];

my %orig2full = ();

open my $FULL, '<', $full_names;
while (my $input = <$FULL>) { 
    chomp $input;
    if ( $input =~ /\A > ([A-Za-z0-9]+ [.]* [_] ([^\s\/]+) ) \/ \d+ [-] \d+ /xms ) { 
        my $full_name = $1;
        my $orig_name = $2;
        $orig2full{$orig_name} = $full_name;
    }
    elsif ( $input =~ /\A > /xms ) { 
        die "Can't parse FASTA header from full-names file $full_names: $input\n";
    }
}
close $FULL;

open my $UPROT, '<', $uprot_dat;
while (my $input = <$UPROT>) {
    chomp $input;
    if ( $input =~ /\A ([A-Z][a-z]+) \t (\S+) \z/xms ) { 
        my $type = $1;
        my $orig_name = $2;
        if ( exists $orig2full{$orig_name} ) {
            print "$type\t$orig2full{$orig_name}\n";
        }
    }
}
close $UPROT;

