#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

my $uni  = q{};
my $tx   = q{};
my $gene = q{};

while (my $input = <>) {
    chomp $input;

    # Sample input:
    # AC   A0A0K0CUW9;
    if ( $input =~ /\A AC \s+ (\S+)[;]\s*\z/xms ) {
         my $new_uni = $1;
         if ($uni) {
             die "Failed to map to ParaSite: $uni\n";
         }
         $uni = $new_uni;
    }

    # Sample input (note that a single $uni can have 2+ of these lines):
    # DR   WBParaSite; ACAC_0000105801-mRNA-1; ACAC_0000105801-mRNA-1; ACAC_0000105801.
    elsif ( $input =~ /\A DR \s+ WBParaSite[;] \s+ (\S+)[;] \s+ \S+[;] \s+ (\S+)[.] \s* \z/xms ) {
        my $new_tx   = $1;
        my $new_gene = $2;

        # error checks
        if (! $uni ) {
            die "Failed to identify UniProt for mapping to \"$new_tx\" and \"$new_gene\"/n";
        }

        # map stuff
        $tx = $new_tx;
        $gene = $new_gene;
        $data_ref->{'uni'}->{$uni}->{'gene'}->{$gene}->{'tx'}->{$tx} = 1;

        # immediately clear $gene and $tx, because a single $uni may have 2+ of these lines
        $tx   = q{};
        $gene = q{};
    }

    # '//' is the end of a normally formatted UniProt record
    elsif ( $input =~ /\A \/ \//xms ) {
        # clear $uni value for the next round
        $uni  = q{};
    }
}

my @unis = sort keys %{ $data_ref->{'uni'} };

foreach my $uni_map (@unis) {
    my @genes = sort keys %{ $data_ref->{'uni'}->{$uni_map}->{'gene'} };
    foreach my $gene_map (@genes) {
        my @txs = sort keys %{ $data_ref->{'uni'}->{$uni_map}->{'gene'}->{$gene_map}->{'tx'} };
        foreach my $tx_map (@txs) {
            print "$uni_map\t$gene_map\t$tx_map\n";
        }
    }    
}

