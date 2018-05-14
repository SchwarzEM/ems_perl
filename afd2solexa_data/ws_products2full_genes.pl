#!/usr/bin/env perl

# ws_products2full_genes.pl -- Erich Schwarz <emsch@caltech.edu>, 11/26/2012.
# Purpose: given WS230-style protein or ncRNA files for C. elegans etc., build a WBGene \t pubname \t seqname table which wbg2fullnames.pl can then use.

use strict;
use warnings;

my $wbgene  = q{};
my $pubname = q{};
my $seqname = q{};

my $data_ref;

while (my $input = <>) {
    if ( $input =~ /\A > /xms ) { 
        if ( $input !~ / \A > \S+ \s .* WBGene\d+ /xms ) { 
            die "Can't parse FASTA header line: $input\n";
        }
        if ( $input =~ / \A > (\S+) \s .* (WBGene\d+) /xms ) {
            $seqname = $1;
            $wbgene  = $2;

            # Deal with two different kinds of CDS suffixes:
            $seqname =~ s/[a-z]\z//;
            $seqname =~ s/\A(\S+\.\d+)\.\d+\z/$1/;

            if ( ( exists $data_ref->{'wbgene'}->{$wbgene}->{'seqname'} ) and ( $seqname ne $data_ref->{'wbgene'}->{$wbgene}->{'seqname'} ) ) { 
                die "Can't decide if seqname for WBGene \"$wbgene\" is old \"$data_ref->{'wbgene'}->{$wbgene}->{$seqname}\" or new \"$seqname\"\n";
            }
            $data_ref->{'wbgene'}->{$wbgene}->{'seqname'} = $seqname;
            $pubname = $seqname;
       }
       if ( $input =~ / locus [:] (\S+) /xms ) { 
           $pubname = $1;
       }
       $data_ref->{'wbgene'}->{$wbgene}->{'pubname'} = $pubname;
   }
}

my @wbgenes = sort keys %{ $data_ref->{'wbgene'} };

foreach my $wbgene1 (@wbgenes) {
    $wbgene  = $wbgene1;
    $pubname = $data_ref->{'wbgene'}->{$wbgene}->{'pubname'};
    $seqname = $data_ref->{'wbgene'}->{$wbgene}->{'seqname'};
    print "$wbgene\t$pubname\t$seqname\n";
}

