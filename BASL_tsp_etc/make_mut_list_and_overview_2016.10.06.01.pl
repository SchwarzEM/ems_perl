#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

my $data_ref;

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

my $overview = "lof_nmd_gene_overview_06oct2016.tsv.txt" ;
$overview = safename($overview);   
open my $OVERVIEW, '>', $overview;

foreach (my $infile = <>) {
    chomp $infile;

    # sample input:
    # lof_nmd_gene_summaries_2016.10.06/jj103.lof_nmd_gene_summary_06oct2016.txt

    if ( $infile =~ /\A \S+ \/ (\S+) \.lof_nmd_gene_summary_06oct2016\.txt \z/xms ) {
        my $genotype = $1;

        my $list = "$genotype.lof_nmd_gene_list_06oct2016.tsv.txt" ;
        $list = safename($list);

        open my $INFILE, '<', $infile;
        while (my $input = <$INFILE>) {
            chomp $input;
            if ( $input =~ / (LOF|NMD) = \( ([^\|\s]+) \| (WBGene\d+) /xms ) {
                my $mut_type  = $1;
                my $gene_name = $2;
                my $wbgene_id = $3;
                my $full_name = "$wbgene_id|$gene_name";
                $data_ref->{'genotype'}->{$genotype}->{'mut_gene'}->{$full_name}->{'mut_type'}->{$mut_type} = 1;
                $data_ref->{'genotype'}->{$genotype}->{'mut_gene'}->{$full_name}->{'wbgene_id'} = $wbgene_id;
                $data_ref->{'genotype'}->{$genotype}->{'mut_type'}->{$mut_type}->{'mut_gene'}->{$full_name} = 1;
            }
        }
        close $INFILE ;
        
        my @mut_genes = sort keys %{ $data_ref->{'genotype'}->{$genotype}->{'mut_gene'} };

        open my $LIST, '>', $list ;
        foreach my $mut_gene (@mut_genes) {
            my $wbgene_id = $data_ref->{'genotype'}->{$genotype}->{'mut_gene'}->{$mut_gene}->{'wbgene_id'} ;

            my $lof_status = q{};
            if ( exists $data_ref->{'genotype'}->{$genotype}->{'mut_type'}->{'LOF'}->{'mut_gene'}->{$mut_gene} ) {
                $lof_status = 'LOF';
            }

            my $nmd_status = q{};
            if ( exists $data_ref->{'genotype'}->{$genotype}->{'mut_type'}->{'NMD'}->{'mut_gene'}->{$mut_gene} ) {
                $nmd_status = 'NMD';
            }

            if ( (! $lof_status ) and (! $nmd_status ) ) {
                die "In genotype $genotype, unable to identify either LOF or NMD for gene $mut_gene\n";
            }

            my $bmp_status = q{};
            if ( exists $known_bmp_genes{$wbgene_id} ) {
                $bmp_status = 'Known_BMP';
            }

            print $LIST "$mut_gene\t$lof_status\t$nmd_status\t$bmp_status\n";
        }
        close $LIST;

        my @lof_genes = ();
        if ( exists $data_ref->{'genotype'}->{$genotype}->{'mut_type'}->{'LOF'}->{'mut_gene'} ) {
            @lof_genes = sort keys %{ $data_ref->{'genotype'}->{$genotype}->{'mut_type'}->{'LOF'}->{'mut_gene'} };
        }

        my @nmd_genes = ();
        if ( exists $data_ref->{'genotype'}->{$genotype}->{'mut_type'}->{'NMD'}->{'mut_gene'} ) { 
            @nmd_genes = sort keys %{ $data_ref->{'genotype'}->{$genotype}->{'mut_type'}->{'NMD'}->{'mut_gene'} };
        }

        if ( (! @lof_genes ) and (! @nmd_genes ) ) {
            die "In genotype $genotype, unable to identify any LOF or NMD for any genes\n";
        }

        my @lof_bmp_genes = ();
        my @nmd_bmp_genes = ();

        if (@lof_genes) {
            foreach my $lof_gene (@lof_genes) {
                my $wbgene_id = $data_ref->{'genotype'}->{$genotype}->{'mut_gene'}->{$lof_gene}->{'wbgene_id'} ;
                if ( exists $known_bmp_genes{$wbgene_id} ) {
                    push @lof_bmp_genes, $lof_gene;
                }
            }
        }
        @lof_bmp_genes = sort @lof_bmp_genes;
        @lof_bmp_genes = uniq(@lof_bmp_genes);

        if (@nmd_genes) {
            foreach my $nmd_gene (@nmd_genes) {
            my $wbgene_id = $data_ref->{'genotype'}->{$genotype}->{'mut_gene'}->{$nmd_gene}->{'wbgene_id'} ;
                if ( exists $known_bmp_genes{$wbgene_id} ) {
                    push @nmd_bmp_genes, $nmd_gene;
                }
            }
        }
        @nmd_bmp_genes = sort @nmd_bmp_genes;
        @nmd_bmp_genes = uniq(@nmd_bmp_genes);

        my $lof_gene_count = @lof_genes;
        my $nmd_gene_count = @nmd_genes;

        my $lof_bmp_gene_count = @lof_bmp_genes;
        my $nmd_bmp_gene_count = @lof_bmp_genes;

        print $OVERVIEW "$genotype\t20000\t$lof_gene_count\t$lof_bmp_gene_count\t$nmd_gene_count\t$nmd_bmp_gene_count\n";
    }
    else { 
        die "Cannot parse input: $infile\n";
    }
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


