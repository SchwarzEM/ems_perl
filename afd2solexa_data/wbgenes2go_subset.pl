#!/usr/bin/env perl

# wbgenes2go_subset.pl -- Erich Schwarz <emsch@its.caltech.edu>, 12/3/2010
# Purpose: given a list of WBGenes and a WormBase gene_association file, generate a list of GO terms linked to them (with or without IEA).

use strict;
use warnings;
use Getopt::Long;

my $wbgene_list = q{};
my $go_assoc    = q{};
my $no_IEA;
my $help;

my $wb_gene = q{};
my $go_term = q{};
my $go_code = q{};

my $data_ref;

GetOptions ( 'wbgene_list=s' => \$wbgene_list,
             'go_assoc=s'    => \$go_assoc,
             'no_IEA'        => \$no_IEA,
             'help'          => \$help,        );

if ($help or (! $wbgene_list) or (! $go_assoc) ) { 
    die "Format: qual2func_table.pl\n",
        "            -w|--wbgene_list [WBGene list, one gene per line]\n",
        "            -g|--go_assoc    [WB GO-association table]\n",
        "            -n|--no_IEA      [ignore IEA annotations]\n",
        "            -h|--help        [print this message]\n",
        "            [print list of GO terms associated with the listed genes to <STDOUT>]\n",
        ;
}

# Sample input lines:
#
# WB      WBGene00000001  aap-1           GO:0008340      WB_REF:WBPaper00005614|PMID:12393910    IMP [...]
# WB      WBGene00000002  aat-1           GO:0006865      PMID:12520011|PMID:12654719     IEA [...]
# Note that we are making no effort to select for annotations belonging only to a single ontology, 
#     e.g., 'P -- but we could if we needed to, by a trivial expansion of the code.

open my $GO_ASSOC, '<', $go_assoc or die "Can't open GO-association table $go_assoc: $!";
while (my $input = <$GO_ASSOC>) { 
    chomp $input;
    if ( $input =~ / \A WB
                     \t (WBGene\d+)
                     \t [^\t]*
                     \t [^\t]*
                     \t (GO[:]\d+)
                     \t [^\t]*
                     \t (I[A-Z]{2})
                     \t /xms        ) {
        $wb_gene = $1;
        $go_term = $2;
        $go_code = $3;
        if ( (! $no_IEA) or ($go_code ne 'IEA') ) { 
            $data_ref->{'gene'}->{$wb_gene}->{'GO'}->{$go_term} = 1;
        }
    }
    else { 
        if ( $input !~ /\A # /xms ) { 
            die "Can't parse input line from $go_assoc: $input\n";
        }
    }
}
close $GO_ASSOC or die "Can't close filehandle to GO-association table $go_assoc: $!";

open my $WBGENE_LIST, '<', $wbgene_list or die "Can't open WBGene ID list $wbgene_list: $!";
while (my $input = <$WBGENE_LIST>) { 
    chomp $input;
    if ( $input =~ /\A (WBGene\d+) /xms ) { 
        $wb_gene = $1;
        if ( exists $data_ref->{'important'}->{$wb_gene} ) { 
            warn "Gene $wb_gene is listed in more than one line of WBGene ID list $wbgene_list!\n";
        }
        $data_ref->{'important'}->{$wb_gene} = 1;
    }
}
close $WBGENE_LIST or die "Can't close filehandle to WBGene ID list $wbgene_list: $!";

foreach my $wb_gene1 (sort keys %{ $data_ref->{'important'} } ) {
    if ( exists $data_ref->{'gene'}->{$wb_gene1}->{'GO'} ) { 
        foreach my $go_term1 ( sort keys %{ $data_ref->{'gene'}->{$wb_gene1}->{'GO'} } ) { 
            print "$go_term1\n";
        }
    }
}

