#!/usr/bin/perl

# ws_ele_prot2tmhmm.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/21/2007.
# Purpose: get elegans proteins/genes for proteins with TMHMM motif.

use strict;
use warnings;
use Ace;

my $database = $ARGV[0];

# Igor Antoschechkin says (10/17/2006):
# AQL gets confused, so make variables unambiguous.

my $aql_query = 
    'select p'                                     # $ref->[0]
    . ', p->Corresponding_CDS'                     # $ref->[1]
    . ', p->Corresponding_CDS->Gene'               # $ref->[2]
    . ', p->Corresponding_CDS->Gene->Public_name'  # $ref->[3]
    . ' from p in class Protein'
    . ' where p->Species like "Caenorhabditis elegans"'
    . ' and p->Feature like "tmhmm"'
   ;

my $db = Ace->connect(-path=>$database) 
    or die "No local database: ",Ace->error;
my @objects = $db->aql($aql_query);

foreach my $ref (sort { $a->[2] cmp $b->[2] } @objects) { 
    if ($ref->[2]) { 
        print "$ref->[2]\t$ref->[3]\t$ref->[1]\t$ref->[0]\n";
    }
}

