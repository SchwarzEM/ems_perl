#!/usr/bin/env perl

# ws_protgene_homols.pl -- Erich Schwarz <emsch@its.caltech.edu>, 4/16/2008.
# Purpose: get named Homol_groups and protein motifs of elegans p-genes.

use strict;
use warnings;
use Ace;

my %results = ();

my $database = $ARGV[0];
if (! $database) { 
    die 'Format: ws_protgene_homols.pl',
        ' [/usr/local/acedb/elegans, or something]',
        "\n",
        ; 
}
unless ( ( -d $database ) and ( -r $database ) ) { 
    die "$database is not readable directory\n"; 
}

# Igor Antoschechkin says (10/17/2006):
# AQL gets confused, so make the variables unambiguous.

my $aql_query = 'select p'                        # $ref->[0]

                . ', h'                           # $ref->[1]
                . ', i'                           # $ref->[2]

                . ', m'                           # $ref->[3]
                . ', o'                           # $ref->[4]

                . ', c'                           # $ref->[5]

                . ', g'                           # $ref->[6]
                . ', n'                           # $ref->[7]
                . ', s'                           # $ref->[8]

                . ' from p in class Protein '

                . ', h in p->Homology_group'
                . ', i in h->Title'

                . ', m in p->Motif_homol'
                . ', o in m->Title'

                . ', c in p->Corresponding_CDS'

                . ', g in c->Gene'
                . ', n in g->Public_name'
                . ', s in g->Sequence_name'

                . ' where p->Species like "Caenorhabditis elegans"'
                ;

my $db = Ace->connect(-path=>$database) 
    or die "No local database: ",Ace->error;

my @aql_refs = $db->aql($aql_query);

my $results_key = '$aql_ref->[0]'
                  . '\t$aql_ref->[1]'
                  . '\t$aql_ref->[2]'
                  . '\t$aql_ref->[3]'
                  . '\t$aql_ref->[4]'
                  . '\t$aql_ref->[5]'
                  . '\t$aql_ref->[6]'
                  . '\t$aql_ref->[7]'
                  . '\t$aql_ref->[8]'
                  ; 

foreach my $aql_ref (@aql_refs) { 
    $results{$results_key} = 1;
}

foreach my $result (sort keys %results) { 
    print "$result\n";
}

