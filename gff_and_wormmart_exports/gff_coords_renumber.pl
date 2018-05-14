#!/usr/bin/perl

# gff_coords_renumber.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/2/2007.
# Purpose: given a nt range and a gene name, extract relative nt coords of CDSes.
# N.B.: output is kludgey and could be improved with better output code; but, it's enough for work now.

use strict;
use warnings;

unless ($#ARGV == 2) { 
    die "Format ./gff_coords_renumber.pl  ",
        "[gff file]  ",
        "[name of sequence to scan]  ",
        "[coords, in desired final orientation]\n",
        ;
}

my $infile    = $ARGV[0];
my $seqname   = $ARGV[1];
my $coords    = $ARGV[2];
my $coord1    = 0;
my $coord2    = 0;
my $index     = 0;
my $direction = 0;

if (! -r $infile) { 
    die "Can't read $infile: $!";
}

if (! $coords =~ /\A \d+ \- \d+ \z /xms) {
    die "Misformatted coordinates: $coords\n";
}

if ($coords =~ /\A (\d+) \- (\d+) /xms) {
    $coord1 = $1;
    $coord2 = $2;

    if ($coord1 == $coord2) { 
        die "$coords has single-residue range\n";
    }

    elsif ($coord1 < $coord2) {
        # Ascending, sense-strand range.
        $direction = 1;
        $index = $coord1 - 1;
    }

    elsif ($coord2 < $coord1) {
        # Descending, antisense-strand range.
        $direction = -1;
        $index = $coord1 + 1;

        # For ease of range-checking later:
        ($coord1, $coord2) = ($coord2, $coord1);
    }

    else { 
        die "Can't make sense of coordinates $coords\n";
    }
}

open my $INFILE, "$infile" or die "Can't open $infile: $!";
while (my $input = <$INFILE>) { 
    chomp $input;
    if ($input =~ / \A            (\S+)       # $sname
                    \s+ .+
                    CDS \s+       (\d+)       # $nt1
                    \s+           (\d+)       # $nt2
                    .+
                    gene_id \s \" ( [^\"]+ )  # $gene
                    \" /xms) { 
        my $sname = $1;
        my $nt1   = $2;
        my $nt2   = $3;
        my $gene  = $4;

        if ( ($sname eq $seqname) 
             and (    ( ($coord1 <= $nt1) and ($nt1 <= $coord2) )
                   or ( ($coord1 <= $nt2) and ($nt2 <= $coord2) ) ) 
           ) {
            $nt1 = redefine_nt($nt1);
            $nt2 = redefine_nt($nt2);
            if ($nt1 > $nt2) { 
                ($nt1, $nt2) = ($nt2, $nt1);
            }
            print "$gene\t$nt1\t$nt2\n";
        }
    }
}
close $INFILE;

sub redefine_nt {
    my $_nt = $_[0];
    if ( ($_nt < $coord1) or ($_nt > $coord2) ) {
        $_nt = 0;
    }
    elsif ( ($coord1 <= $_nt) and ($_nt <= $coord2) ) {
        $_nt = (($_nt - $index) * $direction);
    }
    return $_nt;
}

