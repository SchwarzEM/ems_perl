#!/usr/bin/env perl

# list_aug_gene_complete_status_12jul2016.pl -- Erich Schwarz <ems349@cornell.edu>, 7/12/2016.
# Purpose: given a GFF3-like (not GTF2-like) *.aug output from AUGUSTUS 3.2.2 (or a similar version), tabulate, for each gene, whether it has both a start and a stop codon, only one, or none.
# Note that this is being checked for *transcripts*, but then mapped back to gene names.

# Sample inputs:
#
# nigoni_2016.07.08_001	AUGUSTUS	stop_codon	13553	13555	.	-	0	Parent=nigoni_2016.07.08_001.g1.t1
# nigoni_2016.07.08_001	AUGUSTUS	start_codon	15887	15889	.	-	0	Parent=nigoni_2016.07.08_001.g1.t1
# nigoni_2016.07.08_001	AUGUSTUS	stop_codon	13553	13555	.	-	0	Parent=nigoni_2016.07.08_001.g1.t2
# nigoni_2016.07.08_001	AUGUSTUS	start_codon	15890	15892	.	-	0	Parent=nigoni_2016.07.08_001.g1.t2

use strict;
use warnings;

my $transcript = q{};

my $data_ref;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A \S+ \t AUGUSTUS \t start_codon \t .+ \t Parent= (\S+) \s* \z/xms ) { 
        $transcript = $1;
        $data_ref->{'transcript'}->{$transcript}->{'start_codon'} = 1;
    }
    elsif ( $input =~ /\A \S+ \t AUGUSTUS \t stop_codon \t .+ \t Parent= (\S+) \s* \z/xms ) {
        $transcript = $1;
        $data_ref->{'transcript'}->{$transcript}->{'stop_codon'} = 1;
    }
}

my @transcripts = sort keys %{ $data_ref->{'transcript'} };

foreach my $obs_tx (@transcripts) { 
    if ( $obs_tx =~ /\A (\S+) \.t\d+ \z/xms ) { 
        my $gene1 = $1;
        if ( ( exists $data_ref->{'transcript'}->{$obs_tx}->{'start_codon'} ) and ( exists $data_ref->{'transcript'}->{$obs_tx}->{'stop_codon'} ) ) {
            $data_ref->{'complete_gene'}->{$gene1} = 1;
        }

        if ( (! exists $data_ref->{'transcript'}->{$obs_tx}->{'start_codon'} ) or (! exists $data_ref->{'transcript'}->{$obs_tx}->{'stop_codon'} ) ) {
            $data_ref->{'incomplete_gene'}->{$gene1} = 1;
        }
    } 
    else {
        die "Can't parse transcript name: $obs_tx\n";
    }
}

my @genes = sort keys %{ $data_ref->{'incomplete_gene'} };

foreach my $gene2 (@genes) { 
    if (! exists $data_ref->{'complete_gene'}->{$gene2} ) {
        print "$gene2\n";
    }
}

