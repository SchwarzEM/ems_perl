#!/usr/bin/perl

# largest_wpep_isoform.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/5/2007.
# Purpose:  given a wormpep file (e.g., wormpep170) pick largest isoform per gene (e.g. for InParanoid).
# Usage:    ./largest_wpep_isoform.pl wormpep170 > wormpep170.max.fa

use strict;
use warnings;

my $input       = "";    # input line
my $cds         = "";    # e.g., "T13A10.10a"
my $wprot       = "";    # e.g., "CE30174"
my $gene        = "";    # e.g., "WBGene00000005"
my %prots       = ();    # $cds via "header1", "header2", "wprot" to data [hashref keys]
my %genes2cdses = ();    # all the CDSes for a given gene                 [arrayref keys]

while ($input = <>) { 
    chomp $input;
    if ( ($input =~ /\A >/x) 
             and ($input !~ /\A > .+ status: .+ \z/x) ) { 
        die "Nonparseable header: $input\n";
    }
    if ($input =~ / \A ( > (\S+)     # $2 == $cds
                    \s+ (CE\d+)      # $3 == $wprot
                    \s+ (WBGene\d+)  # $4 == $gene
                    .+ )             # $1 -- insert blurb-stuff after this header-half
                    (status:.+)      # $5 -- and insert it before this one
                    \z /x) { 
        $cds   = $2;
        $wprot = $3;
        $gene  = $4;
        $prots{$cds}->{"header1"} = $1;    # more readable, but *after* '$cds = $2'!
        $prots{$cds}->{"header2"} = $5;
        $prots{$cds}->{"wprot"}  = $wprot; 

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
    # First sort by alphabet, then re-sort by aa_count;
    #     so ties of aa_count pre-resolved.
    my @cdses = sort { $prots{$b}->{"aa_count"} 
                       <=> $prots{$a}->{"aa_count"} } 
                     ( sort @{ $genes2cdses{$gene} } );

    my $key_cds = shift @cdses;

    my @wprots = grep { $_ ne $prots{$key_cds}->{"wprot"} } 
        sort map { $prots{$_}->{"wprot"} } @cdses;

    print $prots{$key_cds}->{"header1"};
    print " [Alt. prots: @wprots] " if (@wprots);
    print $prots{$key_cds}->{"header2"};
    print "\n";

    my $fullseq = $prots{$key_cds}->{"sequence"};
    my @seqs = unpack("a60" x (length($fullseq)/60 + 1), $fullseq );
    foreach my $aa60 (@seqs) {
        print "$aa60\n";
    }
}

