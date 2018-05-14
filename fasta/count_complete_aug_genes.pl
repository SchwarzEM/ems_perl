#!/usr/bin/env perl

# Erich Schwarz <emsch@caltech.edu>, 10/26/2012.
# Purpose: given a GTF-like *.aug output from AUGUSTUS 2.6.1 (or a similar version), count how many genes have both a start and a stop codon.
# Note that this is being checked for *transcripts*, but then mapped back to gene names.

# Sample input:
# Acey_2012.08.05_0001    AUGUSTUS        start_codon     72898   72900   .       +       0       transcript_id "Acey_2012.08.05_0001.g5.t3"; gene_id "Acey_2012.08.05_0001.g5";

use strict;
use warnings;

my $gene       = q{};
my $transcript = q{};

my $data_ref;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A [^\t]* \t AUGUSTUS \t start_codon \t .+ \t [^\t]* transcript_id [ ] \" ([^\"\s]+) \" /xms ) { 
        $transcript = $1;
        $data_ref->{'transcript'}->{$transcript}->{'start_codon'} = 1;
    }
    elsif ( $input =~ /\A [^\t]* \t AUGUSTUS \t stop_codon \t .+ \t [^\t]* transcript_id [ ] \" ([^\"\s]+) \" /xms ) {
        $transcript = $1;
        $data_ref->{'transcript'}->{$transcript}->{'stop_codon'} = 1;
    }
}

my @transcripts = sort keys %{ $data_ref->{'transcript'} };

foreach my $obs_tx (@transcripts) { 
    if ( ( exists $data_ref->{'transcript'}->{$obs_tx}->{'start_codon'} ) and ( exists $data_ref->{'transcript'}->{$obs_tx}->{'stop_codon'} ) ) { 
        if ( $obs_tx =~ /\A (\S+) \.t\d+ \z/xms ) { 
            $gene = $1;
            $data_ref->{'complete_gene'}->{$gene} = 1;
        }
        else { 
            die "Can't parse transcript name: $obs_tx\n";
        }
    }
}

my @genes = sort keys %{ $data_ref->{'complete_gene'} };

foreach my $obs_gene (@genes) { 
    print "$obs_gene\n";
}


