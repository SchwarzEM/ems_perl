#!/usr/bin/perl

# kog_wpep2wpep.pl
# Erich Schwarz <emsch@its.caltech.edu>, 6/6/04

# Purpose: map KOG version of wormpep to modern sequences/names in wormpep125.
# Overall Purpose: allow updating of NCBI KOG_wormpeps so that KOGs can be .ace-ified.

use strict;

unless (($ARGV[1] ne "") && ($ARGV[2] eq "")) { die "Format: kog_wpep2wpep.pl [1line-kog_wpep] [wormpep]\n"; }

my $input_line = "";
my $kog_wpep  = $ARGV[0];
my $wormpep   = $ARGV[1];
my $kog2wpep  = $kog_wpep . ".to." . $wormpep;
my $orphan_kog_wpep = "orphan_" . $kog_wpep . "_vs_" . $wormpep;
my $protname = "";
my $residues = 0;
my %protseq = "";

open WORMPEP, "$wormpep" or die "Can't open $wormpep\n";
while (<WORMPEP>) {
    chomp($input_line = $_);
    if ($input_line =~ /^>(\S+)\s+(CE\d+)\s+locus:(\S+)\s+/) { 
        if ($residues) { $protseq{$protname} = $residues; }
        $protname = "$1\t$2\t$3";
        $residues = "";
    }
    elsif ($input_line =~ /^>(\S+)\s+(CE\d+)\s+/) {
        if ($residues) { $protseq{$protname} = $residues; }
        $protname =  "$1\t$2\t\t";
        $residues = "";
    }
    else {
        $residues = $residues . $input_line;
        $residues =~ s/[\W\d]//g;
    }
}
if ($residues) { $protseq{$protname} = $residues; }
close WORMPEP;

open KOG_WPROTS, "$kog_wpep"               or die "Can't open $kog_wpep\n";
open KOG2WPEP, ">$kog2wpep"                or die "Can't open $kog2wpep\n";
open ORPHAN_KOG_WPEP, ">$orphan_kog_wpep"  or die "Can't open $orphan_kog_wpep\n";
while (<KOG_WPROTS>) {
    my $match_worked = "no";
    chomp ($input_line = $_);
    $input_line =~ m/^(\S+)\s+(\S+)$/ or die "$kog_wpep fails to match expected pattern\n";
    my $kog_protname = $1;
    my $kog_protseq =  $2;
    foreach $protname (sort keys %protseq) {
        if (((index $protseq{$protname}, $kog_protseq) + 1) && ($protseq{$protname}) && ($protname)){ 
            print KOG2WPEP "$kog_protname\t$protname\n";
            $match_worked = "yes";
        }
    }
    if ($match_worked eq "no") {
        print ORPHAN_KOG_WPEP ">$kog_protname\torphan -- found in $kog_wpep but not $wormpep\n";
        while ($kog_protseq =~ /^(\w{1,60})(.*)/)
        {
            my $kog_protseq_printline=$1;
            my $kog_protseq_waiting=$2;
            print ORPHAN_KOG_WPEP "$kog_protseq_printline\n";
            $kog_protseq = $kog_protseq_waiting;
        }
    }
}
close KOG_WPROTS; close KOG2WPEP; close ORPHAN_KOG_WPEP;
