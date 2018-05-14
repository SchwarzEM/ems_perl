#!/usr/bin/perl
# mk_condlist.pl -- Erich Schwarz <emsch@its.caltech.edu>, 10/9/2006.
# Purpose: make a condition/genelist file from a simple gene list.

use strict;
use warnings;

my %cond_genes = ();
my %gene2seq = ();
my %seq2gene = ();
my %gene2cgc = ();
my %cgc2gene = ();

unless ($#ARGV == 2) { die "Format: ./gpep2alitab.pl [gene table] [simple genelist] [condition]\n"; }

open (GENES, $ARGV[0]) or die "Can't open gene table $ARGV[0]: $!";
while (<GENES>) { 
    chomp (my $in1 = $_);
    if ($in1 =~ /(WBGene\d+)\t([^\t]+)\t([^\t]+)\s*$/) { 
        $gene2seq{$1} = $2;
        $seq2gene{$2} = $1;
        $gene2cgc{$1} = $3;
        $cgc2gene{$3} = $1;
    }
    elsif ($in1 =~ /(WBGene\d+)\t([^\t]+)\s*$/) {
        $gene2seq{$1} = $2;
        $seq2gene{$2} = $1;
    }
}
close GENES;

open (GENELIST, $ARGV[1]) or die "Can't open GenPept file $ARGV[0]: $!";
my @genelist = grep {$_ =~ /\S/} <GENELIST>;
chomp @genelist;
my $condition = $ARGV[2];
unless ($condition =~ /\S/) { die "Misformatted condition.\n"; }
print "\"$condition\"\n";
foreach my $id (@genelist) { 
    if ($id =~ /^(WBGene\d+)\D*/) { 
        $cond_genes{$1} = 1;
    }
    elsif (exists $seq2gene{$id}) { 
        $cond_genes{$seq2gene{$id}} = 1;
    }
    elsif (exists $cgc2gene{$id}) {
        $cond_genes{$cgc2gene{$id}} = 1;
    }
    unless ( ($id =~ /^(WBGene\d+)\D*/) or 
             (exists $seq2gene{$id})   or 
             (exists $cgc2gene{$id}) )    { 
        warn "Unrecognized gene name: $id\n"; 
    }
}
close GENELIST;

foreach my $gene (sort keys %cond_genes) { 
    if (exists $gene2cgc{$gene}) { 
        print "$gene" . '|' . "$gene2seq{$gene}" . '|' . "$gene2cgc{$gene}" . "\n";
    }
    elsif (exists $gene2seq{$gene}) {
        print "$gene" . '|' . "$gene2seq{$gene}". "\n";
    }
}

