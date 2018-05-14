#!/usr/bin/env perl

# pfam_hmmscan2annot.pl -- Erich Schwarz <ems394@cornell.edu>, 9/2/2016.
# Purpose: given a gene2cds table and a PFAM/hmmscan 3.0 output, make a 2-column table of genes and PFAM domains.

use strict;
use warnings;
use Getopt::Long;
use autodie;

my $data_ref;

my $gene2cds = q{};
my $pfam     = q{};
my $help;

GetOptions ( 'gene2cds=s' => \$gene2cds,
             'pfam=s'     => \$pfam,
             'help'       => \$help, );

if ( $help or (! $gene2cds) or (! $pfam) ) { 
    die "Format: pfam_hmmscan2annot.pl --gene2cds|-g [gene-to-CDS table] --pfam|-p [PFAM/hmmscan 3.0 tabular output] --help|-h [print this message]\n";
}

open my $GENE, '<', $gene2cds;
while (my $input = <$GENE>) { 
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) { 
        my $gene = $1;
        my $cds  = $2;
        $data_ref->{'cds'}->{$cds}->{'gene'} = $gene;
    }
    else { 
        die "From transcript-to-gene table $gene2cds, can't parse input: $input\n";
    }
}
close $GENE;

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
            $data_ref->{'gene'}->{$gene}->{'pfam_desc'}->{$pfam_desc} = 1;
        }
        else { 
            die "From PFAM table $pfam, can't parse input: $input\n";
        } 
    }
}
close $PFAM;

my @genes = sort keys %{ $data_ref->{'gene'} };

my $header = "Gene\tPFAM\n";

foreach my $gene1 (@genes) {
    my $pfam_text  = q{};
    my @pfam_descs = ();
    if ( exists $data_ref->{'gene'}->{$gene1}->{'pfam_desc'} ) { 
        @pfam_descs = sort keys %{ $data_ref->{'gene'}->{$gene1}->{'pfam_desc'} };
        $pfam_text  = join '; ', @pfam_descs;
    }
    print $header if $header;
    $header = 0;
    print "$gene1\t$pfam_text\n";
}

