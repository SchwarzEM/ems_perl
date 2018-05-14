#!/usr/bin/env perl

# get_wbGO_subset.pl -- Erich Schwarz <emsch@its.caltech.edu>, 7/22/2011.
# Purpose: given a list of "interesting" GO terms and a WormBase GO association file, make a list of genes with any "interesting" terms in a single annotation column.

use strict;
use warnings;

my $relevant_genes = $ARGV[0];
my $interesting    = $ARGV[1];
my $wb_go_assoc    = $ARGV[2];

my $wbgene  = q{};
my $go_term = q{};
my $go_desc = q{};

my $data_ref;

open my $RELEVANT, '<', $relevant_genes or die "Can't open list of relevant genes, $relevant_genes: $!";
while (my $input = <$RELEVANT>) { 
    chomp $input;
    if ( $input =~ /\A (WBGene\d+) /xms ) { 
        $wbgene = $1;
        $data_ref->{'seen_gene'}->{$wbgene} = 1;
    }
}

close $RELEVANT or die "Can't close list of relevant genes, $relevant_genes: $!";

open my $INTERESTING, '<', $interesting or die "Can't open list of interesting GO terms, $interesting: $!";
while (my $input = <$INTERESTING>) { 
    chomp $input;
    # ... protein dephosphorylation	GO:0006470 ... #
    if ( $input =~ / \t ([^\t]+) \t (GO:\d+) \t /xms ) { 
        $go_desc = $1;
        $go_term = $2;
        if ( $go_desc !~ /\S/xms ) { 
            die "Can't parse GO term $go_term, \"$go_desc\"\n";
        }
        $data_ref->{'seen_go_term'}->{$go_term} = $go_desc;
    }
}
close $INTERESTING or die "Can't close filehandle to list of interesting GO terms, $interesting: $!";

open my $WB_GO_ASSOC, '<', $wb_go_assoc or die "Can't open WormBase GO association file, $wb_go_assoc: $!";
while (my $input = <$WB_GO_ASSOC>) {
    chomp $input;
    # WB      WBGene00000001  aap-1           GO:0019901     
    if ( $input =~ /\A WB \t (WBGene\d+) \t [^\t]* \t [^\t]* \t (GO:\d+) \t /xms ) { 
        $wbgene  = $1;
        $go_term = $2;
        if ( ( exists $data_ref->{'seen_go_term'}->{$go_term} ) and ( exists $data_ref->{'seen_gene'}->{$wbgene} ) ) { 
            $data_ref->{'annot_gene'}->{$wbgene}->{'go_term'}->{$go_term} = 1;
        }
    }
}
close $WB_GO_ASSOC or die "Can't close filehandle to WormBase GO association file, $wb_go_assoc: $!";

foreach my $annot_gene ( sort keys %{ $data_ref->{'annot_gene'} } ) { 
    my @go_annots = ();
    foreach my $assoc_go_term ( sort keys %{ $data_ref->{'annot_gene'}->{$annot_gene}->{'go_term'} } ) { 
        my $assoc_go_desc = $data_ref->{'seen_go_term'}->{$assoc_go_term};
        push @go_annots, "$assoc_go_term [$assoc_go_desc]";
    }
    my $go_annot = join '; ', @go_annots;
    print "$annot_gene\t$go_annot\n";
}

