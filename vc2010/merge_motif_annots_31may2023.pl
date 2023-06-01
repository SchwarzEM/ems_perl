#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $stats = q{};
my $mots  = q{};

$stats = $ARGV[0] if $ARGV[0];
$mots  = $ARGV[1] if $ARGV[1];

my $data_ref;

if ( (! $stats ) or (! $mots ) ) {
    die "Format: merge_motif_annots_31may2023.pl [motif stats] [motifs with annots] > [stats with annots] ;\n";
}

if ( (! -r $stats ) or (! -r $mots ) ) {
    die "Cannot read motif stats $stats or motifs with annots $mots\n";
}

open my $STATS, '<', $stats;
while ( my $input = <$STATS> ) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (.+) \z/xms ) {
        my $mot_id    = $1;
        my $mot_stats = $2;
        $data_ref->{'mot_id'}->{$mot_id}->{'mot_stats'} = $mot_stats;
    }
}
close $STATS;

open my $MOTS, '<', $mots;
while ( my $input = <$MOTS> ) {
    # Sample input: MOTIF WBPmat00000001 DAF-16_Furuyama_2000
    chomp $input;
    if ( $input =~ /\A MOTIF \s+ (\S+) \s+ (\S+) \z/xms ) {
        my $mot_id   = $1;
        my $mot_name = $2;
        $data_ref->{'mot_id'}->{$mot_id}->{'mot_name'} = $mot_name;
    }
}
close $MOTS;

my @mot_list = sort keys %{ $data_ref->{'mot_id'} };
foreach my $mot_id (@mot_list) {
    my $mot_stats = q{};
    my $mot_name  = q{};

    # Good to have, but not obligatory:
    if ( exists $data_ref->{'mot_id'}->{$mot_id}->{'mot_name'} ) {
        $mot_name = $data_ref->{'mot_id'}->{$mot_id}->{'mot_name'} ;
    }
    # Put the print command inside this 'if' loop, so that only motifs with statistics get printed:
    if ( exists $data_ref->{'mot_id'}->{$mot_id}->{'mot_stats'} ) {
        $mot_stats = $data_ref->{'mot_id'}->{$mot_id}->{'mot_stats'};
        print "$mot_id\t$mot_name\t$mot_stats\n";
    }
}

