#!/usr/bin/env perl

# extract_GenBank_gbf_prots.pl -- Erich Schwarz <ems@emstech.org>, 10/14/2017.
# Purpose: given a GenBank output from an attempted Sequin loading of an annotated genome, extract the protein sequences into FASTA format (e.g., for cross-checking for accuracy with cd-hit-2b).

use strict;
use warnings;

my $sequence = q{};
my $residues = q{};
my $reading  = 0;

my $data_ref;

while ( my $input = <>) { 
    chomp $input;
    # /product="ACEY_S0270.G837.T1"
    if ( $input =~ /\A \s* \/protein_id=\" \S+ [:] (\S+) \" /xms ) {
        $sequence = $1;
        $residues = q{};
    }

    # deal with rare instances of lines like:  /translation="MPISLPQLPVLVRDEVLQQMNFLE"
    elsif ( $input =~ /\A \s* \/translation=\" (\S+) \" /xms ) {
        $reading = 0;
        my $new_res = $1;
        $residues .= $new_res;

        $data_ref->{'sequence'}->{$sequence}->{'residues'} = $residues;
        $residues = q{};
    }

    # most other lines are like this:  /translation="MPISLPQLPVLVRDEVLQQMNFLE [... etc.]
    elsif ( $input =~ /\A \s* \/translation=\" (\S+)  /xms ) { 
        $reading = 1;
        my $new_res = $1;
        $residues .= $new_res;
    }
    elsif ( $reading and ( $input =~ /\A \s* ([^\"\s]+) \s* \z/xms ) ) { 
        my $new_res = $1;
        $residues .= $new_res; 
    }
    elsif ( $reading and ( $input =~ / \A (.+) \" \s* \z/xms ) ) {
        $reading = 0;
        my $new_res = $1;

        # There may be leading blank residues; indeed; the entire line may be blank until a terminating ".
        $new_res =~ s/\s//g;

        $residues .= $new_res;
        $data_ref->{'sequence'}->{$sequence}->{'residues'} = $residues;
        $residues = q{};
    }
}

my @sequences = sort keys %{ $data_ref->{'sequence'} };
foreach my $sequence1 (@sequences) { 
    my $residues1 = $data_ref->{'sequence'}->{$sequence1}->{'residues'};
    print ">$sequence1\n";
    my @output_lines = unpack("a60" x (length($residues1)/60 + 1), $residues1);
    foreach my $output_line (@output_lines) { 
        if ($output_line =~ /\S/) { 
            print "$output_line\n";
        }
    }
}


