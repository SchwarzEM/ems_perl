#!/usr/bin/perl

# ncbi_idents_from_psiblasts.pl
# Erich Schwarz <emsch@its.caltech.edu>, 10/25/05.
# Purpose: get various idents. from converged psi-blast search of NCBI nr: simple and staved gi and acc. nos., and desc. lines.

use strict;
use warnings;

my $gi_staved_output  = $ARGV[0] . ".gi_staved_nos";
my $gi_simple_output  = $ARGV[0] . ".gi_simple_nos";
my $acc_staved_output = $ARGV[0] . ".acc_staved_nos";
my $acc_simple_output = $ARGV[0] . ".acc_simple_nos";

my $desc_output = $ARGV[0] . ".descrips";

open (GIS_STAVED, ">$gi_staved_output") or die "Can't open GI (staved) number output file: $!";
open (GIS_SIMPLE, ">$gi_simple_output") or die "Can't open GI (simple) number output file: $!"; 
open (ACCS_STAVED, ">$acc_staved_output") or die "Can't open accession (staved) number output file: $!";
open (ACCS_SIMPLE, ">$acc_simple_output") or die "Can't open accession (simple) number output file: $!";

open (DESCS, ">$desc_output") or die "Can't open description output file: $!";

my $reading  = "yes";
my @gi_nos   = ();
my @acc_nos  = ();
my @descrips = ();

while (<>) { 
    chomp(my $input = $_);
    if ($input =~ /^Sequences not found previously or not previously below threshold:/) { 
        $reading = "no";
    }
    elsif ($input =~ /^Sequences producing significant alignments:/) {
        $reading = "yes";
        @gi_nos = ();
        @acc_nos = ();
        @descrips = ();
    }
    elsif ( ( $input =~ /^(gi\|[^\|]+)\|(.+?)[\|]*\s/ ) and ($reading eq "yes") ) {
        push @gi_nos, $1;
        push @acc_nos, $2;
        push @descrips, $input;
    }
    elsif ( ( $input =~ /^(.+?\|[^\|]+)/ ) and ($reading eq "yes") ) {
        push @acc_nos, $1;
    }
    elsif ( $input =~ /CONVERGED/ ) {
        foreach my $gi_no (sort @gi_nos) {
            print GIS_STAVED "$gi_no\n";
            if ($gi_no =~ /^gi\|([^\|]+)[\s\|]*/) { 
                $gi_no = $1;
            }
            print GIS_SIMPLE "$gi_no\n";
        }
        foreach my $acc_no (sort @acc_nos) {
            print ACCS_STAVED "$acc_no\n";
            if ($acc_no =~ /^[^\|]+[\|]+([^\|\.]+)[\.\s\|]*/) { 
                $acc_no = $1;
            }
            print ACCS_SIMPLE "$acc_no\n";
        }
        foreach my $descr (sort @descrips) { 
            print DESCS "$descr\n";
        }
    }
}
