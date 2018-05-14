#!/usr/bin/env perl

# list_internal_stop_cdses.pl -- Erich Schwarz <ems@emstech.org>, 2/22/2014.
# Purpose: given an input stream with CDS DNAs in FASTA format that have unique names, emit a list of those sequences that have an internal stop codon.

use strict;
use warnings;

my $seq = q{};

my %stop = (
    'TAA' => 1,
    'TAG' => 1,
    'TGA' => 1,
);

my $data_ref;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A > (\S+) .* \z/xms ) { 
        $seq = $1;
        if ( exists $data_ref->{'seq'}->{$seq} ) { 
            die "Two different sequences are both, redundantly, called $seq\n";
        }
    }
    elsif ( $input =~ /\A > /xms ) {
        die "Can't parse: $input\n";
    }
    else { 
        $input =~ s/\s//;
        $input =~ tr/[a-z]/[A-Z]/;
        $input =~ tr/U/T/;
        $data_ref->{'seq'}->{$seq}->{'residues'} .= $input;
    }
}

my @seqs = sort keys %{ $data_ref->{'seq'} };
foreach my $seq1 (@seqs) { 
    my $residues = $data_ref->{'seq'}->{$seq1}->{'residues'};

    # For most CDSes, we will have a defined start codon, and can parse them on that basis.
    if ( $residues =~ /\A (ATG (?: \S{3})+ ) (\S{0,2}) \z/xms ) { 
        # If there is a start codon, split into bulk of sequence, which should be in clean triplets, and 0-2 residual 3' nucleotides.
        # Most well-formed CDSes will have zero residual nucleotides, but 
        my $residues_codons   = $1;
        my $residues_residual = $2;

        my $total_length    = length($residues);
        my $residual_length = length($residues_residual);
        if ( $residual_length >= 3 ) { 
            die "Failed to properly parse length of sequence $seq1 (total length $total_length, residual length $residues_residual)\n";
        }

        # Trim off final codon, which will usually (though not always) be a legitimate stop codon.
        if ( $residues_codons =~ /\A (\S{3,}) \S{3} \z/xms ) {
            $residues_codons = $1;
        }
        else { 
            die "For some reason, fail to remove terminal codon of sequence $seq1\n";
        }

        scan_for_stops($residues_codons,$seq1);
    }
    elsif ( $residues =~ /\A (\S{0,2}) ( (?:\S{3})+ ) (\S{3}) \z/xms ) {
        my $residues_residual = $1;
        my $residues_codons   = $2;
        my $final_codon       = $3;

        if ( exists $stop{$final_codon} ) {
            # Enforce correct residue parsing.
            my $total_length    = length($residues);
            my $residual_length = length($residues_residual);
            if ( $residual_length >= 3 ) {
                die "Failed to properly parse length of sequence $seq1 (total length $total_length, residual length $residues_residual)\n";
            }

            scan_for_stops($residues_codons,$seq1);
        }

        else { 
            # This gives a list of those sequences for which there was neither an easy start nor an easy stop codon.
            warn "$seq1\n";
        }
    }
}

my @seqs_w_internal_stops = sort keys %{ $data_ref->{'seq_w_internal_stop'} };
foreach my $seq_w_internal_stop (@seqs_w_internal_stops) { 
    print "$seq_w_internal_stop\n";
}

sub scan_for_stops {
    my $_residues_codons = $_[0];
    my $_seq1            = $_[1];
    while ( $_residues_codons =~ /\A (\S{3}) (.*) \z/xmsg ) {
        my $first_codon  = $1;
        $_residues_codons = $2;
        if ( exists $stop{$first_codon} ) {
            $data_ref->{'seq_w_internal_stop'}->{$_seq1} = 1;
        }
    }
}

