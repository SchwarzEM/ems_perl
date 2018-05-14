#!/usr/bin/env perl

# tabulate_ncoils_x.fa.pl -- Erich Schwarz <ems394@cornell.edu>, 9/28/2014.
# Purpose: given the 'x'-masked NCoils version of a proteome, and its gene-to-CDS mapping, generate a gene-centric table of NCoils annotations.

use strict;
use warnings;
use Getopt::Long;
use autodie;

my $gene2cds = q{};
my $fasta    = q{};

my $gene = q{};
my $cds  = q{};

my $data_ref;

my $help;

GetOptions ( 'gene2cds=s' => \$gene2cds,
             'fasta=s'    => \$fasta, 
             'help'       => \$help, );

if ($help or (! $gene2cds) or (! $fasta) ) { 
    die "Format: tabulate_ncoils_x.fa.pl  --gene2cds|-g  [Gene-to-CDS table] --fasta|-f [ncoils-xed FASTA] --help|-h [print this message]\n";
}

open my $GENE, '<', $gene2cds;
while (my $input = <$GENE>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
        $gene = $1;
        $cds  = $2;
        if ( ( exists $data_ref->{'cds'}->{$cds}->{'gene'} ) and ( $data_ref->{'cds'}->{$cds}->{'gene'} ne $gene ) ) {
            die "CDS $cds is inconsistently mapped to both $data_ref->{'cds'}->{$cds}->{'gene'} and $gene\n";
        }
        $data_ref->{'cds'}->{$cds}->{'gene'} = $gene;
    }
    else {
        die "From gene-to-CDS name table $gene2cds, can't parse input line: $input\n";
    }
}
close $GENE;

open my $FASTA, '<', $fasta;
while (my $input = <$FASTA>) { 
    chomp $input;
    if ( $input =~ /\A > (\S+) /xms ) {
        $cds  = $1;
        $gene = q{};
        if (! exists $data_ref->{'cds'}->{$cds}->{'gene'} ) {
            die "In proteome $fasta, cannot map CDS $cds to gene\n";
        }
        $gene = $data_ref->{'cds'}->{$cds}->{'gene'};
    }
    elsif ( $cds and $gene and ( $input =~ /\S/xms ) ) { 
        $input =~ s/\s//g;
        $data_ref->{'gene'}->{$gene}->{'cds'}->{$cds}->{'seq'} .= $input;
    }
}
close $FASTA or die "Can't close filehandle to  proteome: $fasta\n";

my $header = "Gene\tNCoils\n";

my @genes = sort keys %{ $data_ref->{'gene'} };
foreach my $gene1 (@genes) { 
    my @cdses        = sort keys %{ $data_ref->{'gene'}->{$gene1}->{'cds'} };
    my @coils_annots = ();
    foreach my $cds1 (@cdses) { 
        my $sequence  = $data_ref->{'gene'}->{$gene1}->{'cds'}->{$cds1}->{'seq'};
        my $len_seq   = length($sequence);
        my $len_coils = 0;
        if ( $sequence =~ /[x]+/xms ) { 
            $len_coils = ( $sequence =~ tr/x/x/ );
        }
        my $frac_coils = ($len_coils/$len_seq);
        $frac_coils    = sprintf "%.2f", $frac_coils;
        my $coils_text = "$frac_coils ($len_coils/$len_seq)";
        if ( $frac_coils > 0 ) { 
            push @coils_annots, $coils_text;
        }
    }
    my $coils_annot = join '; ', @coils_annots;
    print $header if $header;
    $header = q{};
    if ( $coils_annot =~ /\S/xms ) { 
        $coils_annot = 'NCoils: ' . $coils_annot;
    }
    print "$gene1\t$coils_annot\n";
}
