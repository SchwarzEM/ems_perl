#!/usr/bin/env perl

use strict ;
use warnings;
use autodie;

my $cds2gene = q{};
my $cds_loj  = q{};

$cds2gene = $ARGV[0] if $ARGV[0];
$cds_loj  = $ARGV[1] if $ARGV[1];

my $data_ref;   # declare as completely unspecified variable, *not* as ' = q{}', or data structures will crash under 'use strict'.

# Format of input $cds2gene: "CDS/transcript name \t gene name"
# Extract this from, e.g., ParaSite protein FASTA header lines such as:
# 
# >augustus-ANCCEYDFT_Contig1-processed-gene-17.0-mRNA-1 transcript=augustus-ANCCEYDFT_Contig1-processed-gene-17.0-mRNA-1 gene=augustus-ANCCEYDFT_Contig1-processed-gene-17.0

if ( (! -e $cds2gene) or (! -e $cds_loj) ) {
    die 'Format: cds_loj_to_genemap_11nov2022.pl [cds2gene] [cds_loj] > [gene1 \t gene2]', "\n";
}

open my $CDS2GENE, '<', $cds2gene;
while (my $input = <$CDS2GENE>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
        my $cds  = $1;
        my $gene = $2;
        $data_ref->{'cds'}->{ $cds }->{'gene'} = $gene;
    }
    else {
        die "From cds2gene $cds2gene, cannot parse input: $input\n";
    }
}
close $CDS2GENE;

# Format of input $cds_loj: excise two columns from a BEDtools loj join using "cut -f 10,20 | grep gene".
# Sample input line:
# ID=cds:maker-ANCCEYDFT_Contig305-pred_gff_snap-mRNA-2.21-mRNA-1;Parent=transcript:maker-ANCCEYDFT_Contig305-pred_gff_snap-mRNA-2.21-mRNA-1;extra_copy_number=0 
# ID=Acey_s0001.v2.g17938.t1.CDS1 gene=Acey_s0001.v2.g17938;Parent=Acey_s0001.v2.g17938.t1 gene=Acey_s0001.v2.g17938;

open my $CDS_LOJ, '<', $cds_loj;
while (my $input = <$CDS_LOJ>) {
    chomp $input;
    if ( $input =~ /\A ID[=]cds[:] ([^\s;]+) .* \t .* gene= ([^\s;]+) .* \z/xms ) {
        my $cds1  = $1;
        my $gene1 = q{};
        my $gene2 = $2;
        if ( exists $data_ref->{'cds'}->{ $cds1 }->{'gene'} ) {
            $gene1 = $data_ref->{'cds'}->{ $cds1 }->{'gene'};
        }
        else {
            # Deal with weird habit of LiftOver: it adds '_\d+' suffixes to some genes it maps.
            $cds1 =~ s/\_\d+\z//;
            if ( exists $data_ref->{'cds'}->{ $cds1 }->{'gene'} ) {
                $gene1 = $data_ref->{'cds'}->{ $cds1 }->{'gene'};
            }
            else {
                warn "Cannot map CDS $cds1 to a gene\n";
                die "Input line: $input\n";
            }
        }
        if ( ( $gene1 =~ /\S/xms ) and ( $gene2 =~ /\S/xms ) ) {
            $data_ref->{'gene1'}->{ $gene1 }->{'gene2'}->{ $gene2 } = 1;
        }
    }
    elsif ( $input !~ /\A .+ \t [.] \z/xms) {
        die "From cds_loj $cds_loj, cannot parse input: $input\n";
    }
}
close $CDS_LOJ;

my @genes1 = grep { /\S/ } sort keys %{ $data_ref->{'gene1'} };
foreach my $gene1 (@genes1) {
    my @genes2 = grep { /\S/ } sort keys %{ $data_ref->{'gene1'}->{ $gene1 }->{'gene2'} };
    foreach my $gene2 (@genes2) {
        print "$gene1\t$gene2\n";
    }
}
