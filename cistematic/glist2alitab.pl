#!/usr/bin/perl
# glist2alitab.pl -- Erich Schwarz <emsch@its.caltech.edu>, 10/2/2006.
# Purpose: read in gene ID list from latest acedb; get gene IDs from condition/genelist file; mark the genes with 1.

use strict;
use warnings;

my %cond_genes = ();
my %gene2seq = ();
my %seq2gene = ();
my %gene2cgc = ();
my %cgc2gene = ();

unless ($#ARGV == 1) { die "Format: ./gpep2alitab.pl [gene table] [condition/genelist file]\n"; }

open (GENES, $ARGV[0]) or die "Can't open gene table $ARGV[0]: $!";
while (<GENES>) { 
    chomp (my $in1 = $_);
    if ($in1 =~ /WBGene\d+\t[^\t]*\-[^\t]*/) { 
        die "Mislocated locus name in gene table: $in1\n";
    }
    elsif ( ($in1 =~ /WBGene\d+\t[^\t]+\t[^\t]+\s*$/) and  # run once, edit when dies
            ($in1 =~ /WBGene\d+\t[^\t]+\t[^\-]+\s*$/) ) { 
        die "Mislocated sequence name in gene table: $in1\n";
    }
    elsif ($in1 =~ /(WBGene\d+)\t([^\t]+)\t([^\t]+)\s*$/) { 
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

open (CONDLIST, $ARGV[1]) or die "Can't open GenPept file $ARGV[0]: $!";
my @condlist = grep {$_ =~ /\S/} <CONDLIST>;
chomp @condlist;
my $condition = shift @condlist;
unless ($condition =~ /^\"[^\"]+\"$/) { die "Misformatted condition line.\n"; }
print "WBGene\t$condition\n";
foreach my $id (@condlist) { 
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
        # This used to be 'die', but that is *far* too fragile to flakes who can't give me accurate names for their genes:
        warn "Unrecognized gene name: $id\n"; 
    }
}
close CONDLIST;

foreach my $gene (sort keys %cond_genes) { 
    if (exists $gene2cgc{$gene}) { 
        print "$gene" . '|' . "$gene2seq{$gene}" . '|' . "$gene2cgc{$gene}" . "\t1\n";
    }
    elsif (exists $gene2seq{$gene}) {
        print "$gene" . '|' . "$gene2seq{$gene}". "\t1\n";
    }
}

