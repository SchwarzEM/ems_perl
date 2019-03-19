#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $gff = q{};
my $ids = q{};
my $bad = q{};

my %good_id = ();

$gff = $ARGV[0] if $ARGV[0];
$ids = $ARGV[1]	if $ARGV[1];
$bad = $ARGV[2] if $ARGV[2];

open my $IDS, '<', $ids;
while ( my $input = <$IDS> ) { 
    chomp $input;
    $good_id{$input} = 1;
}
close $IDS;

open my $GFF, '<', $gff;
while ( my $input = <$GFF> ) {
    chomp $input;
    if ( $input =~ /\A[#]/xms ) {
        print "$input\n";
    }
    elsif ( $input =~ /\A (?:[^\t]* \t){8} gene_id [ ] \" (\S+) \"/xms ) { 
        my $gene_id = $1;
        if ( exists $good_id{$gene_id} and (! $bad ) ) {
            print "$input\n";
        }
        elsif ( $bad and (! exists $good_id{$gene_id} ) ) {
            print "$input\n";
        }
    }
}
close $GFF;
