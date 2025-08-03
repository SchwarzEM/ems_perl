#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $names  = q{};
my $genome = q{};

$names  = $ARGV[0] if $ARGV[0];
$genome = $ARGV[1] if $ARGV[1];

my %rename = ();

if ( (! $names ) or (! $genome ) ) {
    die "Format: rename_gbhbak_03aug2025.pl [renaming table] [genbank genome] > [renamed genome]\n";
}

open my $NAMES, '<', $names;
while ( my $input = <$NAMES> ) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
        my $gb_name = $1;
        my $sg_name = $2;
        if ( exists $rename{$gb_name} ) {
            die "Redundant renaming of $gb_name to both $rename{$gb_name} and $sg_name\n";
        }
        $rename{$gb_name} = $sg_name;
    }
    else {
        die "From renaming table $names, cannot parse: $input\n";
    }
}
close $NAMES;

open my $GENOME, '<', $genome;
while ( my $input = <$GENOME> ) {
    chomp $input;
    if ( $input =~ /\A [>] (\S+) (.*)/xms ) {
        my $gb_name = $1;
        my $comment = $2;
        my $sg_name = q{};
        if ( exists $rename{$gb_name} ) {
            $sg_name = $rename{$gb_name};
            print '>';
            print "$sg_name  $gb_name  $comment\n";
        }
        else {
            print '>';
            print "$gb_name  $comment\n";
        }
    }
    elsif ( $input =~ /\S/xms ) {
        print "$input\n";
    }
}


