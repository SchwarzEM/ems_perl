#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

my $tfseq_to_uniprot = $ARGV[0];
my $tfam_to_tfseq    = $ARGV[1];

open my $SEQ, '<', $tfseq_to_uniprot;
while (my $input = <$SEQ>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
        my $tfseq   = $1;
        my $uniprot = $2;
        $data_ref->{'tfseq'}->{$tfseq}->{'uniprot'}->{$uniprot} = 1;
    }
}
close $SEQ;

open my $TFAM, '<', $tfam_to_tfseq;
while (my $input = <$TFAM>) {
    chomp $input;
    if ( $input =~ /\A (TF\d+) \t (\S+) \z/xms ) {
        my $tfam  = $1;
        my $tfseq = $2;

        # Save time; *only* encode this if there will be an effective mapping later!
        if ( exists $data_ref->{'tfseq'} ) {
            $data_ref->{'tfam'}->{$tfam}->{'tfseq'}->{$tfseq} = 1;
        }
    }
}
close $TFAM;

my @tfams = sort keys %{ $data_ref->{'tfam'} };
foreach my $tfam (@tfams) {
    my @tfseqs = sort keys %{ $data_ref->{'tfam'}->{$tfam}->{'tfseq'} };
    foreach my $tfseq (@tfseqs) {
        my @uniprots = sort keys %{ $data_ref->{'tfseq'}->{$tfseq}->{'uniprot'} };
        foreach my $uniprot (@uniprots) {
            print "$tfam\t$uniprot\n";
        }
    }
}

