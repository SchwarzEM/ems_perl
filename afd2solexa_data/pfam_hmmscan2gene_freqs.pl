#!/usr/bin/env perl

# pfam_hmmscan2annot.pl -- Erich Schwarz <ems394@cornell.edu>, 3/7/2024.
# Purpose: given a cds2gene table and a Pfam/hmmscan 3.0 output, make a 2-column table of Pfam domains and gene frequencies for each domain.

use strict;
use warnings;
use Getopt::Long;
use autodie;

my $data_ref;

my $cds2gene = q{};
my $pfam     = q{};
my $label    = 'Gene_number';
my $help;

GetOptions ( 'cds2gene=s' => \$cds2gene,
             'pfam=s'     => \$pfam,
             'label:s'    => \$label,
             'help'       => \$help, );

if ( $help or (! $cds2gene) or (! $pfam) ) { 
    die "Format: pfam_hmmscan2gene_freqs.pl",
        " --cds2gene|-c [CDS-to-gene table]",
        " --pfam|-p [PFAM/hmmscan 3.0 tabular output]",
        " --label|-l [optional label for gene counts; default is 'Gene_number']",
        " --help|-h [print this message]",
        " > [Pfam domains / gene frequencies, TSV]\n";
}

open my $CDS2GENE, '<', $cds2gene;
while (my $input = <$CDS2GENE>) { 
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) { 
        my $cds  = $1;
        my $gene = $2;
        $data_ref->{'cds'}->{$cds}->{'gene'} = $gene;
    }
    else { 
        die "From transcript-to-gene table $cds2gene, can't parse input: $input\n";
    }
}
close $CDS2GENE;

open my $PFAM, '<', $pfam;
while (my $input = <$PFAM>) {
    chomp $input;
    # Sample input line:
    # Helitron_like_N      PF14214.4  tropicalis_2016.08.11_001.g1.t1 -
    if ( $input !~ /\A [#] /xms ) { 
        if ( $input =~ /\A (\S+) \s+ (\S+) \s+ (\S+) \s+ /xms ) { 
            my $pfam_name = $1;
            my $pfam_acc  = $2;
            my $cds       = $3;
            my $gene      = $data_ref->{'cds'}->{$cds}->{'gene'};
            my $pfam_desc = "$pfam_name [$pfam_acc]";
            $data_ref->{'pfam_desc'}->{$pfam_desc}->{'gene'}->{$gene} = 1
        }
        else { 
            die "From PFAM table $pfam, can't parse input: $input\n";
        } 
    }
}
close $PFAM;

my @pfams = sort keys %{ $data_ref->{'pfam_desc'} };

my $header = "Pfam\t$label\n";

foreach my $indiv_pfam (@pfams) {
    my $gene_count = 0;
    if ( exists $data_ref->{'pfam_desc'}->{$indiv_pfam}->{'gene'} ) {
        my @genes = sort keys %{ $data_ref->{'pfam_desc'}->{$indiv_pfam}->{'gene'} };
        $gene_count = @genes;
    }
    print $header if $header;
    $header = 0;
    print "$indiv_pfam\t$gene_count\n";
}

