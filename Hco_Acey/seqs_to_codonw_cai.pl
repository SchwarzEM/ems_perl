#!/usr/bin/env perl

# seqs_to_codonw_cai.pl -- Erich Schwarz <ems394@cornell.edu>, 1/10/2014.
# Purpose: given a set of CDS DNA sequences, print out a CAI ("weight") value file that can then be used by CodonW to compute CAIs, with user-provided values.

# Notes:
# This is meant to be fed directly into CodonW 1.4.4.
# The file format is obtuse, but can be figured out from the software documentation if one tries hard enough.
# This uses *weights*, which is a mercy, because their computation is considerably simpler than (say) doing 64 weighted Wilcoxons 
#    on codon frequency comparison between highly-expressed genes and all genes in a genome (for Fop, Frequence of optimal codons).
# See indices.txt in CodonW 1.4.4 source code for useful clarifications of CAI.

use strict;
use warnings;

my $sequence        = q{};
my $psuedozero      = 0.01;
my $max_codon_count = 0;

# These are listed in a very particular order that CodonW happens to need.
my @codons = qw(
 TTT TCT TAT TGT 
 TTC TCC TAC TGC 
 TTA TCA TAA TGA 
 TTG TCG TAG TGG 
 CTT CCT CAT CGT 
 CTC CCC CAC CGC 
 CTA CCA CAA CGA 
 CTG CCG CAG CGG 
 ATT ACT AAT AGT 
 ATC ACC AAC AGC 
 ATA ACA AAA AGA 
 ATG ACG AAG AGG 
 GTT GCT GAT GGT 
 GTC GCC GAC GGC 
 GTA GCA GAA GGA 
 GTG GCG GAG GGG 
);

my %cod2fam = (
    "TTT" => "Phe", 
    "TTC" => "Phe", 
    "TTA" => "Leu", 
    "TTG" => "Leu", 
    "CTT" => "Leu", 
    "CTC" => "Leu", 
    "CTA" => "Leu", 
    "CTG" => "Leu", 
    "ATT" => "Ile", 
    "ATC" => "Ile", 
    "ATA" => "Ile", 
    "ATG" => "Met", 
    "GTT" => "Val", 
    "GTC" => "Val", 
    "GTA" => "Val", 
    "GTG" => "Val", 
    "TCT" => "Ser", 
    "TCC" => "Ser", 
    "TCA" => "Ser", 
    "TCG" => "Ser", 
    "CCT" => "Pro", 
    "CCC" => "Pro", 
    "CCA" => "Pro", 
    "CCG" => "Pro", 
    "ACT" => "Thr", 
    "ACC" => "Thr", 
    "ACA" => "Thr", 
    "ACG" => "Thr", 
    "GCT" => "Ala", 
    "GCC" => "Ala", 
    "GCA" => "Ala", 
    "GCG" => "Ala", 
    "TAT" => "Tyr", 
    "TAC" => "Tyr", 
    "TAA" => "STOP", 
    "TAG" => "STOP", 
    "CAT" => "His", 
    "CAC" => "His", 
    "CAA" => "Gln", 
    "CAG" => "Gln", 
    "AAT" => "Asn", 
    "AAC" => "Asn", 
    "AAA" => "Lys", 
    "AAG" => "Lys", 
    "GAT" => "Asp", 
    "GAC" => "Asp", 
    "GAA" => "Glu", 
    "GAG" => "Glu", 
    "TGT" => "Cys", 
    "TGC" => "Cys", 
    "TGA" => "STOP", 
    "TGG" => "Trp", 
    "CGT" => "Arg", 
    "CGC" => "Arg", 
    "CGA" => "Arg", 
    "CGG" => "Arg", 
    "AGT" => "Ser", 
    "AGC" => "Ser", 
    "AGA" => "Arg", 
    "AGG" => "Arg", 
    "GGT" => "Gly", 
    "GGC" => "Gly", 
    "GGA" => "Gly", 
    "GGG" => "Gly", 
);

my $data_ref;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A > \S+ /xms ) { 
        if ($sequence) { 
            process_sequence($sequence);
        }
        $sequence = q{};
    }
    else { 
        $input =~ s/\s//g;
        $sequence .= $input;
    }
}
process_sequence($sequence);

foreach my $codon (@codons) {
    my $family = $cod2fam{$codon};
    my $indiv_codon_count = $data_ref->{'codon'}->{$codon}->{'count'};
    if (      (! exists $data_ref->{'family'}->{$family}->{'max_codon_count'}               ) 
          or ( $indiv_codon_count > $data_ref->{'family'}->{$family}->{'max_codon_count'} ) ) { 
        $data_ref->{'family'}->{$family}->{'max_codon_count'} = $indiv_codon_count;
    }
}

foreach my $codon (@codons) {
    my $family = $cod2fam{$codon};
    my $indiv_codon_count = $data_ref->{'codon'}->{$codon}->{'count'};
    $max_codon_count      = $data_ref->{'family'}->{$family}->{'max_codon_count'};
    my $codon_weight = ( $indiv_codon_count / $max_codon_count );
    if ( $codon_weight < $psuedozero ) {
        $codon_weight = $psuedozero;
    }

    # All codon weights should look like '0.2677518'.
    $codon_weight = sprintf("%.7f", $codon_weight);

    print "$codon_weight\n";
}

sub process_sequence { 
    my $_sequence = $_[0];

    # Make all nucleic acid sequences into ACGT.
    $_sequence =~ tr/uU/tT/;
    $_sequence =~ tr/acgt/ACGT/;

    # Reject non-NA sequences.
    if ( $_sequence =~ /[^ACGT]/xms ) { 
        die "Can't parse sequence $_sequence\n";
    }

    while ( $_sequence =~ /\A([ACGT]{3})(.*)\z/xmsg ) {
        my $_codon  = $1;
        $_sequence  = $2;
        $data_ref->{'codon'}->{$_codon}->{'count'}++;
    }
    return;
}
