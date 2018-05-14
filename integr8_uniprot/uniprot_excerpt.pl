#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;

# Program: uniprot_excerpt.pl -- Erich Schwarz <emsch@its.caltech.edu>, 1/23/2006
# Purpose: Get subsets of *.dat files from Integr8; also export a single, sorted list of taxonomy etc. for that subset.

my %protlist = ();

my $protfile = basename(shift @ARGV);   # first argument should be list of proteins; all others can be *.dat
my $prot_ex  = $protfile . ".dat";
my $prot_spp = $protfile . ".spp";
open (PROTFILE, "$protfile")   or die "Can't open protein list $protfile: $!";
open (PROT_EX,  ">$prot_ex")  or die "Can't open protein .dat excerpt $prot_ex : $!";
open (PROT_SPP, ">$prot_spp") or die "Can't open protein taxonomy list $prot_spp: $!";
while (<PROTFILE>) { 
    chomp (my $input = $_);
    $protlist{$input} = 1;
}
close PROTFILE;

my $do_excerpt = 0;
my $sourcefile = "";
my $protid     = "";
my $species    = "no_species_name";
my $taxon      = "no_taxon_name";
my $ncbitax    = "";
my %proteomes  = ();

while (<>) {
    chomp (my $input = $_);
    if ($do_excerpt) { print PROT_EX "$input\n"; }
    $sourcefile = basename($ARGV);
    if ($input =~ /^ID\s+(\S+)/) { 
        $protid = $1;
        if (exists $protlist{$protid}) { 
            print PROT_EX "\n$input\n";
            $do_excerpt = 1;
        }
        $taxon   = "no_taxon_name";
        $species = "no_species_name";
        $ncbitax = "";
    }
    elsif ($input =~ /^OC\s+(.+)\.\s*$/) {
        if ($taxon eq "no_taxon_name")     { $taxon  = ""; }
        unless ($taxon eq "no_taxon_name") { $taxon .= $1; }
    }
    elsif ($input =~ /^OC\s+(.+[^.])\s*$/) {
        if ($taxon eq "no_taxon_name")     { $taxon = ""; }
        unless ($taxon eq "no_taxon_name") { $taxon = ($taxon . $1 . ' '); }
    }
    elsif ($input =~ /^OS\s+(.+)\.\s*$/) {
        if ($species eq "no_species_name")     { $species  = ""; }
        unless ($species eq "no_species_name") { $species .= $1; }
    }
    elsif ($input =~ /^OS\s+(.+[^.])\s*$/) {
        if ($species eq "no_species_name")     { $species = ""; } 
        unless ($species eq "no_species_name") { $species = ($species . $1 . ' '); }
    }
    elsif ($input =~ /^OX\s+(NCBI_TaxID=\d+)/) { 
        $ncbitax = $1;
    }
    elsif ($input =~/^SQ\s+/) { 
        if ($do_excerpt) { $proteomes{$sourcefile} = "$taxon\t$species\t\t$ncbitax"; }
    }
    elsif ($input =~ /^\/\//) { 
        if ($do_excerpt) { print PROT_EX "\n"; }
        $do_excerpt = 0;
    }
}
close PROT_EX;

foreach $sourcefile (sort { $proteomes{$a} cmp $proteomes{$b} } keys %proteomes) { 
    print PROT_SPP "$proteomes{$sourcefile}\t$sourcefile\n";
}
close PROT_SPP;

