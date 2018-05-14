#!/usr/bin/perl

# ws_gene_names_refs.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/16/2007.
# Purpose: get briggsae protein-to-gene table from command-line-specified local ACeDB database.

use strict;
use warnings;
use Ace;

my $database = $ARGV[0];

my $aql_query = "select p"                    # $ref->[0]
                . ", p->Corresponding_CDS"    # $ref->[1]
                . " from p in class Protein "
                . " where p->Species like \"Caenorhabditis briggsae\""
   ;

my $db = Ace->connect(-path=>$database) or die "No local database: ",Ace->error;
my @objects = $db->aql($aql_query);
foreach my $ref (sort { $a->[0] cmp $b->[0] } @objects) { 
    print "$ref->[0]\t$ref->[1]\n";
}

