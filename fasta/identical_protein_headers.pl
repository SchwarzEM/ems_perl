#!/usr/bin/env perl

# ident_proteome_members.pl -- Erich Schwarz <ems394@cornell.edu>, 11/5/2018.
# Purpose: given two proteomes, print TSV with headers of first proteome's members followed by all headers of matching second proteome members.
# Note that this treats the first proteome as an index and can assign 2+ hits in the second proteome to each first-proteome member.

use strict;
use warnings;
use autodie;

my $data_ref;

my $first  = q{};
my $second = q{};

$first  = $ARGV[0] if $ARGV[0];
$second = $ARGV[1] if $ARGV[1];

my $seq  = q{};
my $head = q{};
my $res  = q{};

open my $FIRST, '<', $first;
while (my $input = <$FIRST>) {
    chomp $input;
    if ( $input =~ /\A [>] ((\S+) .*) \z/xms ) {
        my $new_head = $1;
        my $new_seq  = $2;
        my $new_res  = q{};
        if ( $seq and $head and $res ) {
            $data_ref->{'residues'}->{$res}->{'first_seq'}->{$seq}->{'header'} = $head;
        }
        $seq  = $new_seq;
        $head = $new_head;
        $res  = $new_res;
    }
    elsif ( $input =~ /\A [>] /xms ) {
        die "Cannot parse input: $input\n";
    }
    elsif ( $input =~ /\S/xms ) { 
        $input =~ s/\s//g; 
        $res = $res . $input;
    }
}
if ( $seq and $head and $res ) {
    $data_ref->{'residues'}->{$res}->{'first_seq'}->{$seq}->{'header'} = $head;
}
close $FIRST;

$seq  = q{};
$head = q{};
$res  = q{};

open my $SECOND, '<', $second;
while (my $input = <$SECOND>) {
    chomp $input;
    if ( $input =~ /\A [>] ((\S+) .*) \z/xms ) {
        my $new_head = $1;
        my $new_seq  = $2;
        my $new_res  = q{};
        if ( $seq and $head and $res ) {
            $data_ref->{'residues'}->{$res}->{'second_seq'}->{$seq}->{'header'} = $head;
        }
        $seq = $new_seq;
        $head = $new_head;
        $res = $new_res;
    }
    elsif ( $input =~ /\A [>] /xms ) {
        die "Cannot parse input: $input\n";
    }
    elsif ( $input =~ /\S/xms ) { 
        $input =~ s/\s//g; 
        $res = $res . $input;
    }
}
if ( $seq and $head and $res ) {
    $data_ref->{'residues'}->{$res}->{'second_seq'}->{$seq}->{'header'} = $head;
}
close $SECOND;

my @raw_seqs = sort keys %{ $data_ref->{'residues'} };
foreach my $raw_seq (@raw_seqs) {
    if ( ( $data_ref->{'residues'}->{$raw_seq}->{'first_seq'} ) and ( exists $data_ref->{'residues'}->{$raw_seq}->{'second_seq'} ) ) {
        my @first_seqs = sort keys %{ $data_ref->{'residues'}->{$raw_seq}->{'first_seq'} };
        foreach my $first_seq (@first_seqs) {
            my $first_header = $data_ref->{'residues'}->{$raw_seq}->{'first_seq'}->{$first_seq}->{'header'};
            print "$first_header";

            my @second_seqs = sort keys %{ $data_ref->{'residues'}->{$raw_seq}->{'second_seq'} };
            foreach my $second_seq (@second_seqs) {
                my $second_header = $data_ref->{'residues'}->{$raw_seq}->{'second_seq'}->{$second_seq}->{'header'};            
                print "\t$second_header";
            }
            print "\n";
        }
    }
}

