#!/usr/bin/perl

# proteome_psisubset.pl:  Erich Schwarz <emsch@its.caltech.edu>, 4/4/06.
# Purpose: from psi-blast scan of specialized proteome, extract subset of hits and add header tags (e.g., species data).

use strict;
use warnings;
use File::Basename;

unless (($#ARGV == 2) and (-f $ARGV[0]) and (-f $ARGV[1])) { die "Format: hack1.pl [psi-blast] [proteome] \"quoted tag\" \n"; }

my $blastfile = $ARGV[0];
my $proteome  = $ARGV[1];
my $headtag   = $ARGV[2];
my $warnings  = basename($proteome) . ".warnings";
my $reading  = "no";
my %psi_names = ();

open (BLAST, "$blastfile") or die "Can't open psi-blast file $blastfile: $!";
while (<BLAST>) { 
    chomp(my $input = $_);
    if ($input =~ /^Sequences not found previously or not previously below threshold:/) { 
        $reading = "no";
    }
    elsif ($input =~ /^Sequences used in model and found again:/) { 
# was                /^Sequences producing significant alignments:/
        $reading = "yes";
        %psi_names = ();
    }
    elsif ( ( $input =~ /^(\S+)\s*.+/ ) and ($reading eq "yes") ) { 
        $psi_names{$1} = 1;
    }
    elsif ( $input =~ /CONVERGED/ ) { 
        $reading = "no";
        last; 
    }  # leave %psi_names populated
}
close BLAST;

open (PROT, "$proteome") or die "Can't open proteome file $proteome: $!";
while (<PROT>) {
    chomp (my $fasta_line = $_);
    if ( ($fasta_line =~ /^>(\S+)\s*/) and (! exists $psi_names{$1}) ) {
        $reading = "no";
    }
    elsif ( ($fasta_line =~ /^>(\S+)\s*/) and (exists $psi_names{$1}) and ($psi_names{$1} == 1) ) {
        print "$fasta_line  $headtag\n";
        $reading = "yes";
        $psi_names{$1} = 0;
    }
    elsif ($reading eq "yes") {
        print "$fasta_line\n";
    }
}
$reading = "no";
close PROT ;

my @warns = sort {$b <=> $a} values %psi_names;
if ($warns[0]) { 
    open (WARNINGS, ">$warnings") || die "Can't open $warnings. $!\n";
    foreach my $key (sort keys %psi_names) {
        if ($psi_names{$key} == 1) { 
            print WARNINGS "$key not found in $proteome\n";
        }
    }
}
