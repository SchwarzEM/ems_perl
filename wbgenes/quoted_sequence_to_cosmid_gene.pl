#!/usr/bin/perl

use strict;

# Program: quoted_sequence_to_cosmid_gene.pl
# Erich Schwarz <emsch@its.caltech.edu>, 10/21/02
# 
# Purpose: Get individual cosmid.number genes from things like "C07A4.7a".

my $infile = "";

if ($#ARGV != 0) 
{
    print "Input file? ";
    chomp ($infile = <STDIN>);
} 
else 
{
    $infile = $ARGV[0];
}

my $rough_outfile = ($infile . ".rough_outfile");
my $outfile       = ($infile . ".cosmid_genes");

open (INFILE, "$infile") || die "Wormpep file $infile not found. $!\n";
open (ROUGH_OUTFILE, ">$rough_outfile") || die "Couldn't open file $rough_outfile. $!\n";

while (<INFILE>) 
{
    if ($_ =~ /^\"(\S+)[a-z]{1,}\"/)
    {
        print ROUGH_OUTFILE "$1\n";
    }
    elsif ($_ =~ /^\"(\S+)\"/) 
    {
        print ROUGH_OUTFILE "$1\n";
    }
}

close INFILE;
close ROUGH_OUTFILE;

system "sort $rough_outfile | uniq - > $outfile";
