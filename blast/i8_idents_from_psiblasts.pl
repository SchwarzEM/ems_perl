#!/usr/bin/perl

# i8_idents_from_psiblasts.pl
# Erich Schwarz <emsch@its.caltech.edu>, 1/20/06.
# Purpose: get either single or triple id. nos. from converged psi-blast search of uniprot_to_clean_tfa.pl output.

use strict;
use warnings;
use File::Basename;

my $filename = basename($ARGV[0]);

my $i8_names  = $filename . ".i8_names";
my $i8_4names = $filename . ".i8_4names";

open (I8_NAMES,  ">$i8_names")  or die "Can't open integr8 names list: $!";
open (I8_4NAMES, ">$i8_4names") or die "Can't open integr8 four-names list: $!";

my $reading  = "yes";
my @acc_nos  = ();

while (<>) { 
    chomp(my $input = $_);
    if ($input =~ /^Sequences not found previously or not previously below threshold:/) { 
        $reading = "no";
    }
    elsif ($input =~ /^Sequences producing significant alignments:/) {
        $reading = "yes";
        @acc_nos = ();
    }
    elsif ( ( $input =~ /^(\S+(\s+\S+){3})/ ) and ($reading eq "yes") ) {
        push @acc_nos, $1;
    }
    elsif ( $input =~ /CONVERGED/ ) {
        foreach my $acc_no (sort @acc_nos) {
            $acc_no =~ m/^(\S+)\s+/;
            print I8_NAMES "$1\n";
            print I8_4NAMES "$acc_no\n";
        }
    }
}

