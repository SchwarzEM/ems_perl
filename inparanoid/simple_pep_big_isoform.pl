#!/usr/bin/perl

# simple_wpep_big_isoform.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/6/2007.
# Purpose:  given a (worm|brig|rem)pep file, and species name, pick largest isoform per gene (e.g. for InParanoid).
# Usage:    ./simple_pep_big_isoform.pl brigpep2g briggsae > brigpep2.max.fa

use strict;
use warnings;

unless ($#ARGV == 1) {
    die "Format: ./simple_pep_big_isoform.pl",
        " [proteome FASTA] [species]",
        " > [size-selected proteome]\n";
}

my $input       = "";
my $cds         = "";
my $gene        = "";
my $species     = pop @ARGV;
my %prots       = ();
my %genes2cdses = ();

my %ok_spp = ( elegans  => 1,
               briggsae => 1,
               remanei  => 1, );
if (! $ok_spp{$species} ) { 
    die "Species name \"$species\" not recognized\n";
}

while ($input = <>) { 
    chomp $input;
    if ( ($species eq 'elegans')
          and ($input =~ /\A (\S+) .+ (WBGene\d+)/x) ) {
        ($cds, $gene) = ($1, $2);
        $prots{$cds}->{"header"} = $input;
        push @{ $genes2cdses{$gene} }, $cds;
    }
    if ( ($species eq 'briggsae') 
          and ($input =~ /\A > (CBP\d+) \s+ (CBG\d+) /x) ) { 
        ($cds, $gene) = ($1, $2);
        $prots{$cds}->{"header"} = $input;
        push @{ $genes2cdses{$gene} }, $cds;
    }
    if ( ($species eq 'remanei') 
          and ($input =~ /\A > ( (\S+) \. \d+ ) \s+ /x) ) {
        ($cds, $gene) = ($1, $2);
        $prots{$cds}->{"header"} = $input;
        push @{ $genes2cdses{$gene} }, $cds;
    }
    if ( ($input !~ /\A >/) and ($input =~ /[a-zA-Z]/) ) { 
        $input    =~ tr/[^a-zA-Z]//;
        $input    =~ tr/[a-z]/[A-Z]/;
        my $aa_no = ($input =~ tr/[A-Z]/[A-Z]/); # tr/x/x/ as counting.
        $prots{$cds}->{"sequence"} .= $input;
        $prots{$cds}->{"aa_count"} += $aa_no;
    }
}

foreach $gene (sort keys %genes2cdses ) { 
    my @cdses = sort { $prots{$b}->{"aa_count"} <=> $prots{$a}->{"aa_count"} } 
                ( sort @{ $genes2cdses{$gene} } );
    # First sort by alphabet, then re-sort by aa_count; so ties of aa_count pre-resolved.
    print "DEBUG: ", $cdses[0], " -- ", $prots{$cdses[0]}->{"header"}, "\n";
    my $seq = $prots{$cdses[0]}->{"sequence"};
    my @seqs = unpack("a60" x (length($seq)/60 + 1), $seq );
    foreach $seq (@seqs) {
        print "$seq\n";
    }
}

