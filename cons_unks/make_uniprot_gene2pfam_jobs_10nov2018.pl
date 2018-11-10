#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my @inputs1 = qw(
    ../uniprot_idmappings/arabidopsis_uniprot2gene_2018.11.05.tsv.txt
    ../uniprot_idmappings/dictyostelium_uniprot2gene_2018.11.05.tsv.txt
    ../uniprot_idmappings/drosophila_uniprot2gene_2018.11.05.tsv.txt
    ../uniprot_idmappings/elegans_uniprot2gene_2018.11.05.tsv.txt
    ../uniprot_idmappings/human_uniprot2gene_2018.11.05.tsv.txt
    ../uniprot_idmappings/mouse_uniprot2gene_2018.11.05.tsv.txt
    ../uniprot_idmappings/pombe_uniprot2gene_2018.11.05.tsv.txt
    ../uniprot_idmappings/pristionchus_uniprot2gene_2018.11.05.tsv.txt
    ../uniprot_idmappings/yeast_uniprot2gene_2018.11.05.tsv.txt
    ../uniprot_idmappings/zebrafish_uniprot2gene_2018.11.05.tsv.txt 
);

my @inputs2 = qw(
    ../uniprot_idmappings/arabidopsis_uniprot2pfam.tsv.txt
    ../uniprot_idmappings/dictyostelium_uniprot2pfam.tsv.txt
    ../uniprot_idmappings/drosophila_uniprot2pfam.tsv.txt
    ../uniprot_idmappings/elegans_uniprot2pfam.tsv.txt
    ../uniprot_idmappings/human_uniprot2pfam.tsv.txt
    ../uniprot_idmappings/mouse_uniprot2pfam.tsv.txt
    ../uniprot_idmappings/pombe_uniprot2pfam.tsv.txt
    ../uniprot_idmappings/pristionchus_uniprot2pfam.tsv.txt
    ../uniprot_idmappings/yeast_uniprot2pfam.tsv.txt
    ../uniprot_idmappings/zebrafish_uniprot2pfam.tsv.txt 
);

my $input3 = '../pfam_idmappings/pfam_acc2name_10nov2018.txt';

my @taxa = qw(
    arabidopsis
    dictyostelium
    drosophila
    elegans
    human
    mouse
    pombe
    pristionchus
    yeast
    zebrafish 
);

my $max_val = @inputs1;
$max_val--;
my @vals = (0..$max_val);

foreach my $i (@vals) {
    print '    $SCRATCH/ems_perl/cons_unks/uniprot2gene.family_09nov2018.pl ';
    print "$inputs1[$i] $inputs2[$i] $input3";
    print ' 1>';
    print "$taxa[$i]"; 
    print "_gene2pfam_10nov2018.tsv.txt";
    print ' 2>';
    print "$taxa[$i]"; 
    print "_gene2pfam_10nov2018.err ;\n";
}


