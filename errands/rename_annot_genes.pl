#!/usr/bin/perl

# rename_annot_genes.pl
# Erich Schwarz <emsch@its.caltech.edu>, 7/23/04 
# Purpose: rename Locus, CDS, and Sequence to Gene in hum.-readable annot.s, with geneIDs.WSxxx

use strict;

unless ($#ARGV == 1) { die "Format: ./rename_annot_genes.pl [geneIDs.WSxxx] [annots]\n"; } 
my $genelist = $ARGV[0]; 
my $annot = $ARGV[1]; 
my $new_annot = $annot . ".renamed"; 
my %gene_idents = ""; 
my $input = ""; 
my $gene = "";
my $oldname = ""; 
my $printing = 0;

open GENE_LIST, "$genelist" or die "Can't open $genelist: $!"; 
foreach (<GENE_LIST>) {
    chomp ($input = $_);
    $input =~ s/,/\t/g;
    while ($input =~ /^(\S+)(\s.*?)(\S+)\s*$/) {
        $gene_idents{$3} = $1;
        $input = $1 . $2;
    }
} 
close GENE_LIST;

# open DEBUG, ">debug_list" or die "Can't open debugging list: $!";
# foreach my $key (sort keys %gene_idents) {
#     print DEBUG "$key\t$gene_idents{$key}\n";
# }

open ANNOT_FILE, "$annot"  or die "Can't open $annot: $!"; 
open NEW_ANNOT, ">$new_annot" or die "Can't open $new_annot: $!";

while (<ANNOT_FILE>) {
    chomp($input = $_);
    if ($input =~ /^(CDS|Sequence)\s+:\+\"\S+[b-z]{1}\"/) {
        $printing = 0;
    }
    elsif ($input =~ /^(Locus|CDS|Sequence)\s+:\s+\"(\S+)\"/) { 
        $oldname = $2;
        $oldname =~ s/a$//;
        if ((!$gene_idents{$oldname}) && (!$oldname =~ /[b-z]{1}$/)) {
            warn "Failure to recognize old name $oldname; translated as \"$gene_idents{$oldname}\"?\n";
            $printing = 0;
        }
        elsif ($gene_idents{$oldname}) {
            $printing = 1;
            print NEW_ANNOT "\n";
            print NEW_ANNOT 'Gene : "';
            print NEW_ANNOT $gene_idents{$oldname};
            print NEW_ANNOT '"   // ';
            print NEW_ANNOT "$oldname\n";
        }
    }
    elsif ($input =~ /^(\w+_description\s\")(\S+)\"/) {
        my $part1 = $1;
        $oldname = $2;
        my $part3 = $';
        if ($printing == 1) { 
            $oldname =~ s/a$//;
            if ($gene_idents{$oldname}) { print NEW_ANNOT "$part1" . "$gene_idents{$oldname}" . "\"$part3\n"; }
        }
    }
    elsif ($input =~ /^LongText : \"(\S+)\"/) {
        $oldname = $1;
        if ($printing == 1) {
            $oldname =~ s/a$//;
            if ($gene_idents{$oldname}) { print NEW_ANNOT "LongText : \"$gene_idents{$oldname}\"\n"; }
        }
    }
    elsif ($input =~ /^\*{3}LongTextEnd\*{3}/) { 
        if ($printing == 1) { print NEW_ANNOT "$input\n\n"; }
        $printing = 0;
    }
    elsif ($printing == 1) { print NEW_ANNOT "$input\n"; } 
} 
close ANNOT_FILE; close NEW_ANNOT;
