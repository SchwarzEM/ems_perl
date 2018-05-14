#!/usr/bin/env perl

# orth_remgene2protlist.pl -- Erich Schwarz <emsch@its.caltech.edu>, 4/22/2008.
# Purpose: from gene-centric OrthoMCLs and rempep, get protein list.

use strict;
use warnings;
use Carp;

unless ($#ARGV == 1) { 
    die "Format: ./orth_remgene2protlist.pl [gene-centric OrthoMCL list] [rempep]\n";
}

my $ortho_list = $ARGV[0];
my $rempep     = $ARGV[1];
my %remgenes   = ();
my %remprots   = ();

open my $ORTHO_LIST, '<', $ortho_list
  or croak "Can't open gene-centric OrthoMCL file $ortho_list: $!";

while ( my $input = <$ORTHO_LIST> ) {
    if ( $input =~ / (\S+) \(remanei\) /xms ) {
        my $gene = $1;
        $remgenes{$gene} = 1;
    }
}
close $ORTHO_LIST;

open my $REMPEP, '<', $rempep
  or croak "Can't open remanei proteome $rempep: $!";

while ( my $input = <$REMPEP> ) {
    if ( $input =~ /\A > ( (\S+) \.\d+ ) \s /xms ) {
        my $protein = $1;
        my $gene    = $2;
        if ( exists $remgenes{$gene} ) {
            $remprots{$protein} = 1;
        }
    }
}
close $REMPEP;

foreach my $remprot ( sort keys %remprots ) {
    print "$remprot\n";
}

