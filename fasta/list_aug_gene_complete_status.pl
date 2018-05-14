#!/usr/bin/env perl

# list_aug_gene_complete_status.pl -- Erich Schwarz <emsch@caltech.edu>, 3/4/2013.
# Purpose: given a GTF-like *.aug output from AUGUSTUS 2.6.1 (or a similar version), tabulate, for each gene, whether it has both a start and a stop codon, only one, or none.
# Note that this is being checked for *transcripts*, but then mapped back to gene names.

# Sample input:
# Acey_2012.08.05_0001    AUGUSTUS        start_codon     72898   72900   .       +       0       transcript_id "Acey_2012.08.05_0001.g5.t3"; gene_id "Acey_2012.08.05_0001.g5";

use strict;
use warnings;

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
    if ( $obs_tx =~ /\A (\S+) \.t\d+ \z/xms ) { 
        my $gene1 = $1;
        if ( ( exists $data_ref->{'transcript'}->{$obs_tx}->{'start_codon'} ) and ( exists $data_ref->{'transcript'}->{$obs_tx}->{'stop_codon'} ) ) {
            $data_ref->{'complete_gene'}->{$gene1}->{'status'}->{'Complete'} = 1;
        }
        elsif ( exists $data_ref->{'transcript'}->{$obs_tx}->{'start_codon'} ) {
            $data_ref->{'complete_gene'}->{$gene1}->{'status'}->{'Start codon'} = 1;
        }
        elsif ( exists $data_ref->{'transcript'}->{$obs_tx}->{'stop_codon'} ) {
            $data_ref->{'complete_gene'}->{$gene1}->{'status'}->{'Stop codon'} = 1;
        }
        else {
            $data_ref->{'complete_gene'}->{$gene1}->{'status'}->{'No start/stop codons'} = 1;
        }
    } 
    else {
        die "Can't parse transcript name: $obs_tx\n";
    }
}

my $header = "Gene\tStart/stop codons\n";

my @genes = sort keys %{ $data_ref->{'complete_gene'} };

foreach my $gene2 (@genes) { 
    # Print header once, if there's any data to head:
    print $header if $header;
    $header = q{};
    my $status = q{};

    # For each gene, figure out how best to summarize it.
    if ( exists $data_ref->{'complete_gene'}->{$gene2}->{'status'}->{'Complete'} ) {
        $status = 'Complete';
    }
    elsif ( ( exists $data_ref->{'complete_gene'}->{$gene2}->{'status'}->{'Start codon'} ) or ( exists $data_ref->{'complete_gene'}->{$gene2}->{'status'}->{'Stop codon'} ) ) { 
        my @traits = ();
        if ( exists $data_ref->{'complete_gene'}->{$gene2}->{'status'}->{'Start codon'} ) {
            push @traits, 'Start codon';
        }
        if ( exists $data_ref->{'complete_gene'}->{$gene2}->{'status'}->{'Stop codon'} ) {
            push @traits, 'Stop codon';
        }
        $status = join '; ', @traits;
    }
    else { 
        if (! $data_ref->{'complete_gene'}->{$gene2}->{'status'}->{'No start/stop codons'} ) { 
            die "Failed to make sense of whether gene $gene2 was incomplete or not!\n";
        }
        $status = 'No start/stop codons';
    }
    print "$gene2\t$status\n";
}


