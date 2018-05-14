#!/usr/bin/env perl

# ws_gene2exp.pl -- Erich Schwarz <emsch@its.caltech.edu>, 4/18/2008.
# Purpose: get genes and expression patterns for elegans.

use strict;
use warnings;
use Ace;
use Carp;

my %wbgenes_seen = ();
my %wbgene2clust = ();
my %wbgene2stage = ();
my %wbgene2expre = ();

my $database = $ARGV[0];
if ( !$database ) {
    croak 'Format: ws_gene2exp.pl',
          ' [/usr/local/acedb/elegans, or something]',
          "\n",
          ;
}
if ( (! -d $database ) || (! -r $database ) ) {
    croak "$database is not readable directory\n";
}

# Igor Antoschechkin says (10/17/2006):
# AQL gets confused, so make the variables unambiguous.

my $aql_query = 'select g'    # $ref->[0] -- GeneID
  . ', c'                     # $ref->[1] -- Expression_cluster
  . ', e'                     # $ref->[2] -- Expr_pattern (passed over)
  . ', s'                     # $ref->[3] -- Life_stage
  . ', a'                     # $ref->[4] -- Anatomy_term
  . ', t'                     # $ref->[5] -- Term

  . ' from g in class Gene '           # g! not p! in class..., and ->Species!
  . ', c in g->Expression_cluster'
  . ', e in g->Expr_pattern'
  . ', s in e->Life_stage'
  . ', a in e->Anatomy_term'
  . ', t in a->Term'

  . ' where g->Species like "Caenorhabditis elegans"';

my $db = Ace->connect( -path => $database )
  or croak 'No local database: ', Ace->error;

my @aql_refs = $db->aql($aql_query);

foreach my $aql_ref (@aql_refs) {
    my $wbgene     = $aql_ref->[0];
    my $expr_clust = $aql_ref->[1];
    my $life_stage = $aql_ref->[2];

    my $anat_term = $aql_ref->[4];
    my $anat_desc = $aql_ref->[5];
    $anat_term = "$anat_term [$anat_desc]";

    $wbgenes_seen{$wbgene} = 1;
    if ( $expr_clust =~ /\S/xms ) {
        $wbgene2clust{$wbgene}->{$expr_clust} = 1;
    }
    if ( $life_stage =~ /\S/xms ) {
        $wbgene2stage{$wbgene}->{$life_stage} = 1;
    }
    if ( $anat_term =~ /\S/xms ) {
        $wbgene2expre{$wbgene}->{$anat_term} = 1;
    }
}

foreach my $wbgene ( sort keys %wbgenes_seen ) {
    print $wbgene;

    print "\t";
    if ( keys %{ $wbgene2clust{$wbgene} } ) {
        my $clusters = sort_join_hrefkeys( $wbgene2clust{$wbgene} );
        print $clusters;
    }

    print "\t";
    if ( keys %{ $wbgene2stage{$wbgene} } ) {
        my $life_stages = sort_join_hrefkeys( $wbgene2stage{$wbgene} );
        print $life_stages;
    }

    print "\t";
    if ( keys %{ $wbgene2expre{$wbgene} } ) {
        my $anat_terms = sort_join_hrefkeys( $wbgene2expre{$wbgene} );
        print $anat_terms;
    }
}

sub sort_join_hrefkeys {
    my $hash_ref = $_[0];
    my @terms    = sort keys %{$hash_ref};
    my $termline = join q{; }, @terms;
    return $termline;
}

