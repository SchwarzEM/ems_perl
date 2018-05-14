#!/usr/bin/perl

# count_annot_generefs.pl; Erich Schwarz, <emsch@its.caltech.edu>, 3/26/04.
# Purpose: tabulate quantity and quality of annotations from tabular output of sort_annot_generefs.pl
# Usage: count_annot_generefs.pl <[input_file(s) or <STDIN>] >[output_file or <STDOUT>] ; can be used in pipeline.
# Three lines of example input (would be a very tiny, tab-delimited input):

# 20030124	annots-22jan2003_old-format.ace	abf-6	pmid8940016	01/24/2003
# 20030915	annots-15sep2003_old-format.ace	abf-6	pmid8940016	09/15/2003
# 20040317	annots-17mar2004_old-format.ace	abf-6	pmid8940016	03/17/2004

use strict;

my $input_line     = "";

my $date           = "unset";
my $gene           = "";
my $ref            = "";
my $filename       = "";
my $cal_date       = "";
my $gene_assoc     = "";

my %gene_count     = "";
my %ref_count      = "";
my %annot_count    = "";

my %old_genes      = "";
my %new_genes      = "";
my %reannot_genes  = "";

my %old_assocs     = "";
my %new_assocs     = "";

my $gen_num        = "";
my $ref_num        = "";
my $ann_num        = "";
my $reann_num      = 0;
my $reann_status   = "start";
 
foreach $input_line (<>) {
    chomp ($input_line);
    my @annot_data = split /\t/, $input_line;   # $date, $filename,  $gene, $ref, $cal_date
    unless ($filename == $annot_data[1] && $date == $annot_data[0]) { 
        unless ($date eq "unset") {
            $gen_num   = keys %gene_count;
            $ref_num   = keys %ref_count;
            $ann_num   = $annot_count{$filename};
            if ($reann_status eq "unset") {           # after entering very first annots
                %reannot_genes = "";
                $reann_status = "set";
                $reann_num = 0;
            }
            if ($reann_status eq "start") {           # before entering first annots
                %reannot_genes = "";
                $reann_status  = "unset";
                $reann_num = 0;
            } 
            unless (($reann_num eq "start") || ($reann_num eq "unset")) {
                $reann_num = keys %reannot_genes;  # finally, reannotations can be counted
            }
            print "$cal_date\t$gen_num\t$ref_num\t$ann_num\t$reann_num\t$filename\n";
            my @new_old_genes = keys %new_genes;
            foreach my $new_gene (@new_old_genes) {
                $old_genes{$new_gene}++;
            }
            %new_genes = "";
            @new_old_genes = "";
            my @new_old_assocs = keys %new_assocs;
            foreach my $new_assoc (@new_old_assocs) {
                $old_assocs{$new_assoc}++;
            }
            %new_assocs = "";
            @new_old_assocs = "";
        }
        %gene_count    = "";
        %ref_count     = "";
        %annot_count   = "";
        %new_genes     = "";
        %new_assocs    = "";
        $filename      = "";
        # do NOT reset %reannot_genes, %old_genes, or %old assocs !
    }
    $date        = $annot_data[0];
    $filename    = $annot_data[1];
    $gene        = $annot_data[2];
    $ref         = $annot_data[3];
    $cal_date    = $annot_data[4];
    $gene_assoc  = $annot_data[2] . "\t" . $annot_data[3];
    unless (exists $old_genes{$gene}) { 
        $new_genes{$gene}++;
    }
    if (exists $old_genes{$gene}) {
        unless (exists $old_assocs{$gene_assoc}) { 
            unless ($reann_num eq "unset") {
                $reannot_genes{$gene}++;
                $new_assocs{$gene_assoc}++;
            }
            if ($reann_num eq "unset") {
                $old_assocs{$gene_assoc}++;
            }
        }
    }
    $annot_count{$filename}++;
    $gene_count{$gene}++;
    $ref_count{$ref}++;
}

$gen_num   = keys %gene_count;
$ref_num   = keys %ref_count;
$ann_num   = $annot_count{$filename};
$reann_num = keys %reannot_genes;
print "$cal_date\t$gen_num\t$ref_num\t$ann_num\t$reann_num\t$filename\n";  # clear last entry
