#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $key_fasta  = $ARGV[0];
my $full_table = $ARGV[1];

my $data_ref;

open my $FASTA, '<', $key_fasta;
while (my $input = <$FASTA>) {
    chomp $input;
    if ( $input =~ /\A > ([^\s\/]+) \/ \d+ [-] \d+ /xms ) { 
        my $seq = $1;
        $data_ref->{'seen'}->{$seq} = 1;
    }
    elsif ( $input =~ /\A > /xms ) {
        die "From key FASTA $key_fasta, can't parse header: $input\n";
    }
}
close $FASTA;

open my $FULL, '<', $full_table;
while (my $input = <$FULL>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t .* \z/xms ) { 
        my $seq = $1;
        if ( exists $data_ref->{'read'}->{$seq} ) {
            die "Redundant data entry for sequence $seq: $input\n";
        }
        if ( exists $data_ref->{'seen'}->{$seq} ) {
            $data_ref->{'seq'}->{$seq}->{'data'} = $input;
            $data_ref->{'read'}->{$seq} = 1;
            delete $data_ref->{'seen'}->{$seq};
        }
    }
}
close $FULL;

my @data_seqs = sort keys %{ $data_ref->{'seq'} };
foreach my $seq (@data_seqs) {
    my $data = $data_ref->{'seq'}->{$seq}->{'data'};
    print "$data\n";
}

my @unannotated_seqs = sort keys %{ $data_ref->{'seen'} };
foreach my $unannot_seq (@unannotated_seqs) {
    warn "Failed to find annotation data for: $unannot_seq\n";
}

