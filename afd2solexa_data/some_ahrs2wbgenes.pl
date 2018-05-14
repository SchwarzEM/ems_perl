#!/usr/bin/env perl

# some_ahrs2wbgenes.pl -- Erich Schwarz <emsch@its.caltech.edu>, 5/26/2011.
# Purpose: given a set of genes screened with Ahringer IDs, Ahringer's mapping, and a canonical gene list, get a set of genes with full gene names.

use strict;
use warnings;

my $list  = $ARGV[0];
my $map   = $ARGV[1];
my $names = $ARGV[2];

my ( $cds, $cgc, $cloneid, $wbgene);

my %name2gene  = ();
my %clone2gene = ();
my %seen       = ();

# Sample input:
# WBGene00006385	R119.6	taf-4

open my $NAMES, '<', $names or die "Can't open name table $names: $!";
while (my $input = <$NAMES>) {
    chomp $input; 
    if ( $input =~ / \A (WBGene\d+) \t ([^\t]+) \t ([^\t]+) /xms ) { 
        $wbgene = $1;
        $cds    = $2;
        $cgc    = $3;
        my $wbgene_orig = $wbgene;
        $wbgene = $wbgene . '|' . $cds . '|' . $cgc;
        $name2gene{$wbgene_orig} = $wbgene;
    }
    elsif ( $input =~ / \A (WBGene\d+) \t ([^\t]+) \t /xms ) {
        $wbgene = $1;
        $cds    = $2;
        my $wbgene_orig = $wbgene;
        $wbgene = $wbgene . '|' . $cds ;
        $name2gene{$wbgene_orig} = $wbgene;
    }
    else { 
        die "Can't parse input line from name table $names: $input!\n";
    }
}
close $NAMES or die "Can't close filehandle to name table $names: $!";

# Sample input:
# I-1E17	R119.6	sjj_R119.6	WBGene00006385	R119.6	CDS

open my $MAP, '<', $map or die "Can't open Ahringer clone/gene name map $map: $!";
while (my $input = <$MAP>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \s .+ \s (WBGene\d+) \s /xms ) { 
        $cloneid = $1;
        $wbgene  = $2;
        if (! exists $name2gene{$wbgene} ) { 
            die "Can't map gene name of: $input\n";
        }
        $wbgene = $name2gene{$wbgene};
        $clone2gene{$cloneid} = $wbgene;
    }
    else { 
        die "Can't parse, from $map, input line: $input\n";
    }
}
close $MAP or die "Can't close filehandle to Ahringer clone/gene name map $map: $!";

# Sample input:
# taf-4	I-1E17

open my $LIST, '<', $list or die "Can't open list of selected gene/Ahringer clones $list: $!";
while (my $input = <$LIST>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \s+ (\S+) /xms ) { 
        $cloneid = $2;
        if ( exists $clone2gene{$cloneid} ) {
            print "$clone2gene{$cloneid}\n";
        }
        else { 
            warn "Can't parse, from $list, input line: $input\n";
        }
    }
}
close $LIST or die "Can't close filehandle to list of selected gene/Ahringer clones $list: $!";

