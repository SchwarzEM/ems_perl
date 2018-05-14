#!/usr/bin/perl

# ws_gene_names.pl -- Erich Schwarz <emsch@its.caltech.edu>, 12/14/2006.
# Purpose: get gene / mol. name / locus name list.

use strict;
use warnings;
use Ace;

my %results = ();

my $database = $ARGV[0];
unless ($database) { 
    die "Format: ws_gene_names.pl [/usr/local/acedb/elegans, or something]\n"; 
}
unless ((-d $database) and (-r $database)) { 
    die "$database is not readable directory\n"; 
}

my $aql_query = "select g"                  # $aql_ref->[0]
                . ", g->Sequence_name"      # $aql_ref->[1]
                . ", g->CGC_name"           # $aql_ref->[2]
                . " from g in class Gene "
                . " where g->Species like \"Caenorhabditis elegans\""
   ;

my $db = Ace->connect(-path=>$database) or die "No local database: ",Ace->error;
my @aql_refs = $db->aql($aql_query);
foreach my $aql_ref (@aql_refs) { 
    $results{"$aql_ref->[0]\t$aql_ref->[1]\t$aql_ref->[2]"} = 1;
}
foreach my $result (sort keys %results) { 
    print "$result\n";
}

