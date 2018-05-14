#!/usr/bin/env perl

# merge_wbgene_lists.pl -- Erich Schwarz <emsch@its.caltech.edu>, 4/18/2008.
# Purpose: WBGene-st. list 1 (sorted) + list 2  -> list 1, but w/ 2's full info replacing WBGene names.

use strict;
use warnings;
use Carp;

my $list_1 = $ARGV[0];
my $list_2 = $ARGV[1];

my %wbgene2info = ();

open my $ADDED_INFO, '<', $list_2
  or croak "Can't open file $list_2: $!";
while ( my $input = <$ADDED_INFO> ) {
    chomp $input;
    if ( $input =~ / \A ((WBGene\d+).+)  /xms ) {
        my $info   = $1;
        my $wbgene = $2;
        $wbgene2info{$wbgene} = $info;
    }
}
close $ADDED_INFO;

open my $ORIG_LIST, '<', $list_1
  or croak "Can't open file $list_1: $!";
while ( my $input = <$ORIG_LIST> ) {
    if ( $input =~ / \A ( (WBGene\d+) \S+ ) (.*?) \z  /xms ) {
        my $orig_wbgene_id = $1;
        my $wbgene         = $2;
        my $orig_info      = $3;
        if ( exists $wbgene2info{$wbgene} ) {
            print $wbgene2info{$wbgene};
        }
        if ( !exists $wbgene2info{$wbgene} ) {
            print $orig_wbgene_id;
        }
        print $orig_info;
    }
}
close $ORIG_LIST;

