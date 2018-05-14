#!/usr/bin/perl

# largest_pep_isoform.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/13/2007.
# Purpose:  given a (worm|brig|rem)pep file, and species name, pick largest isoform/gene (e.g. for InParanoid).
# Usage:    ./largest_pep_isoform.pl wormpep170 elegans > wormpep170.max.fa

use strict;
use warnings;

unless ($#ARGV == 1) {
    die "Format: ./largest_pep_isoform.pl",
        " [proteome FASTA] [species]",
        " > [size-selected proteome]\n";
}

my $input       = "";         # input line
my $cds         = "";         # e.g., "T13A10.10a"
my $wprot       = "";         # e.g., "CE30174"
my $gene        = "";         # e.g., "WBGene00000005"
my $species     = pop @ARGV;  # e.g., "elegans"
my %prots       = ();         # $cds via "header1", "header2", "wprot" to data [hashref keys]
my %genes2cdses = ();         # all the CDSes for a given gene                 [arrayref keys]

my %ok_spp = ( elegans  => 1,
               briggsae => 1,
               remanei  => 1, );
if (! $ok_spp{$species} ) {
    die "Species name ", $species, " not recognized\n";
}

while ($input = <>) { 
    chomp $input;

    if (   ($species eq "elegans")
       and ($input =~ / \A ( > (\S+)     # $2 == $cds
                        \s+ (CE\d+)      # $3 == $wprot
                        \s+ (WBGene\d+)  # $4 == $gene
                        .+ )             # $1 -- insert blurb-stuff after this header-half
                        (status:.+)      # $5 -- and insert it before this one
                        \z /x) ) {       # /x regex
        $cds   = $2;
        $wprot = $3;
        $gene  = $4;
        $prots{$cds}->{"header1"} = $1;    # more readable, but *after* '$cds = $2'!
        $prots{$cds}->{"header2"} = $5;
        $prots{$cds}->{"wprot"}  = $wprot; 

        push @{ $genes2cdses{$gene} }, $cds;
    }

    elsif ( ($species eq 'briggsae')
          and ($input =~ /\A > (CBP\d+)     # $1 == $wprot
                          \s+  (WBGene\d+)  # $2 == $gene
                          /x) ) {
        $wprot = $1;
        $gene  = $2;
        $cds   = $wprot;
        $prots{$cds}->{"header1"} = $input;
        $prots{$cds}->{"wprot"}   = $wprot;

        push @{ $genes2cdses{$gene} }, $cds;
    }

    elsif ( ($species eq 'remanei')
          and ($input =~ /\A > ((\S+)      # $2 == $gene
                               \.\d+)      # $1 == $wprot
                               \s+ /x) ) {
        $wprot = $1;
        $gene  = $2;
        $cds   = $wprot;
        $prots{$cds}->{"header1"} = $input;
        $prots{$cds}->{"wprot"}  = $wprot;

        push @{ $genes2cdses{$gene} }, $cds;
    }

    elsif ( ($input !~ /\A >/) and ($input =~ /[a-zA-Z]/) ) {
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

    # Don't bother if no alt. isoforms left!
    my @wprots = ();
    if (@cdses) { 
        @wprots = grep { $_ ne $prots{$key_cds}->{"wprot"} } 
                  sort map { $prots{$_}->{"wprot"} } @cdses;
    }

    print $prots{$key_cds}->{"header1"};
    print " [Alt. prots: @wprots] " if (@wprots);
    print $prots{$key_cds}->{"header2"} 
        if ($prots{$key_cds}->{"header2"}); 
    print "\n";

    my $fullseq = $prots{$key_cds}->{"sequence"};
    my @seqs = unpack("a60" x (length($fullseq)/60 + 1), $fullseq );
    foreach my $aa60 (@seqs) {
        print "$aa60\n" if ($aa60 =~ /\S/);
    }
}

