#!/usr/bin/perl

# largest_wpep_isoform.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/1/2007.
# Purpose:  given a wormpep file (e.g., wormpep170) pick largest isoform per gene (e.g. for InParanoid).
# Usage:    ./largest_wpep_isoform.pl wormpep170 > wormpep170.max.fa

use strict;
use warnings;

my $input       = "";
my $cds         = "";
my $gene        = "";
my %prots       = ();
my %genes2cdses = ();

while ($input = <>) { 
    chomp $input;
    if ($input =~ /^>(\S+).+(WBGene\d+)/) { 
        ($cds, $gene) = ($1, $2);
        $prots{$cds}->{"header"} = $input;
        push @{ $genes2cdses{$gene} }, $cds;
    }
    elsif ($input =~ /[a-zA-Z]/) { 
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
    print $prots{$cdses[0]}->{"header"}, "\n";
    my $seq = $prots{$cdses[0]}->{"sequence"};
    my @seqs = unpack("a60" x (length($seq)/60 + 1), $seq );
    foreach $seq (@seqs) {
        print "$seq\n";
    }
}

