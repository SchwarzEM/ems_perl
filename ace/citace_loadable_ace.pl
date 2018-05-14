#!/usr/bin/perl

# citace_loadable_ace.pl
# Purpose: make .ace files load-proof for citace (touch object, blank object, then remake object).
# Erich Schwarz <emsch@its.caltech.edu>, 10/6/04

use strict;

my $recording = 0;
my @record = ();
my %wipetag = ();
$^I = ".txt";

while (<>) {
    my $input = $_;
    if ($input =~ /^\S+\s+:\s+\"\S+\"/) { 
        $recording = 1;
        push @record, $input;
    }
    elsif (($input =~ /\S+/) and ($recording)) {
        push @record, $input;
        if ($input =~ /^(\S+)\s*/) {
            $wipetag{$1} = 1;
        }
    }
    elsif ($recording) {
        $recording = 0;
        print "\n";
        print $record[0];
        foreach my $wipekey (sort keys %wipetag) {
            print "-D ". $wipekey . "\n";
        }
        print "\n";
        foreach (@record) { print; }
        print "\n";
        @record = ();
    }
}
        
