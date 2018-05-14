#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

my $header = "Mut_gene\tBMP_status\tGenotype_count\tGenotypes_affected\n";

my %known_bmp_genes = (
    WBGene00000900 => 1,
    WBGene00000936 => 1,
    WBGene00003055 => 1,
    WBGene00003056 => 1,
    WBGene00004856 => 1,
    WBGene00004857 => 1,
    WBGene00004858 => 1,
    WBGene00004860 => 1,
    WBGene00004862 => 1,
    WBGene00006324 => 1,
    WBGene00006776 => 1,
    WBGene00011437 => 1,
    WBGene00020649 => 1,
    WBGene00022154 => 1,
    WBGene00023491 => 1,
);

my $overview = "lof_genes_06oct2016.tsv.txt" ;
$overview = safename($overview);   

while (my $infile = <>) {
    chomp $infile;

    # sample input:
    # lof_nmd_gene_summaries_2016.10.06/jj103.lof_nmd_gene_summary_06oct2016.txt

    if ( $infile =~ /\A \S+ \/ (\S+) \.lof_nmd_gene_summary_06oct2016\.txt \z/xms ) {
        my $genotype = $1;

        open my $INFILE, '<', $infile;
        while (my $input = <$INFILE>) {
            chomp $input;
            if ( $input =~ / LOF = \( ([^\|\s]+) \| (WBGene\d+) /xms ) {
                my $gene_name = $1;
                my $wbgene_id = $2;
                my $full_name = "$wbgene_id|$gene_name";
                $data_ref->{'mut_gene'}->{$full_name}->{'genotype'}->{$genotype} = 1;
                $data_ref->{'mut_gene'}->{$full_name}->{'wbgene_id'} = $wbgene_id;
            }
        }
        close $INFILE ;
    }
    else {
        die "Cannot parse input file: $infile\n";
    }
}

my @mut_genes = sort keys %{ $data_ref->{'mut_gene'} };

open my $OVERVIEW, '>', $overview;
foreach my $mut_gene (@mut_genes) {
    my $wbgene_id = $data_ref->{'mut_gene'}->{$mut_gene}->{'wbgene_id'};

    my $bmp_status = q{};
    if ( exists $known_bmp_genes{$wbgene_id} ) {
        $bmp_status = 'Known_BMP';
    }

    my @genotypes = sort keys %{ $data_ref->{'mut_gene'}->{$mut_gene}->{'genotype'} };
    my $genotype_count = @genotypes;
    my $genotype_list  = join ', ', @genotypes;

    print $OVERVIEW $header if $header;
    $header = q{};

    print $OVERVIEW "$mut_gene\t$bmp_status\t$genotype_count\t$genotype_list\n";
}
close $OVERVIEW;

sub safename {
    my $filename = $_[0];
    my $orig_filename = $filename;
    if (-e $orig_filename) {
        my $suffix1 = 1;
        $filename = $filename . ".$suffix1";
        while (-e $filename) {
            $suffix1++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix1";
        }
    }
    return $filename;
}


