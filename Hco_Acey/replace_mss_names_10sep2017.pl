#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $mss_table = q{};
my $proteome  = q{};

my %tx2mss = ();

$mss_table = $ARGV[0];
$proteome  = $ARGV[1];

open my $MSS_TABLE, '<', $mss_table;
while (my $input = <$MSS_TABLE>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) /xms ) { 
        my $tx  = $1;
        my $mss = $2;
        if ( exists $tx2mss{$tx} ) {
            die "Transcript $tx is mapped to MSSes $mss and $tx2mss{$tx}\n";
        }
        $tx2mss{$tx} = $mss;
    }
    else {
        die "In MSS table $mss_table, cannot parse: $input\n";
    }
}
close $MSS_TABLE;

open my $PROTEOME, '<', $proteome;
while (my $input = <$PROTEOME>) {
    chomp $input;
    if ( $input =~ / [>] (\S+) (\b.*) \z/xms ) { 
        my $tx    = $1;
        my $annot = $2;
        if ( exists $tx2mss{$tx} ) {
            my $mss = $tx2mss{$tx};
            print '>', "$mss  $tx  $annot\n";
        }
        else {
            print "$input\n";
        }
    }
    else {
        print "$input\n";
    }
}
close $PROTEOME;

