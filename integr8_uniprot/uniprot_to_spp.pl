#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;

# Program: uniprot_to_clean_tfa.pl -- Erich Schwarz <emsch@its.caltech.edu>, 1/20/2006
# Purpose: Convert *.dat files from Integr8 into a single, sorted list of taxonomy, species, NCBI ID no., and sourcefile.

my $sourcefile = "";
my $species    = "no_species_name";
my $taxon      = "no_taxon_name";
my $ncbitax    = "";
my %proteomes  = ();

while (<>) {
    $sourcefile = basename($ARGV);
    chomp (my $input = $_);
    if ($input =~ /^ID\s+(\S+)/) { 
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
        $proteomes{$sourcefile} = "$taxon\t$species\t\t$ncbitax";
    }
}

foreach $sourcefile (sort { $proteomes{$a} cmp $proteomes{$b} } keys %proteomes) { 
    print "$proteomes{$sourcefile}\t$sourcefile\n";
}

