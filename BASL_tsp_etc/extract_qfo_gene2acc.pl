#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $ok_ids   = $ARGV[0];
my $gene2acc = $ARGV[1];

my $data_ref;

open my $OK_IDS, '<', $ok_ids;
while (my $input = <$OK_IDS>) {
    chomp $input;
    if ( $input =~ /\A \S+ \z/xms ) {
        $data_ref->{'ok_id'}->{$input} = 1;
    }
    else {
        die "From OK IDs file $ok_ids, can't parse: $input\n";
    }
}
close $OK_IDS;

open my $GENE2ACC, '<', $gene2acc;
while (my $input = <$GENE2ACC>) {
    chomp $input;
    if ( $input =~ /\A [^\t]+ \t (\S+) \z/xms ) {
        my $uniprot_id = $1;
        if ( exists $data_ref->{'ok_id'}->{$uniprot_id} ) {
            print "$input\n";
        }
    }
    else {
        die "From gene2acc file $gene2acc, can't parse: $input\n";
    }
}
close $GENE2ACC;

