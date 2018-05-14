#!/usr/bin/env perl

# qual2func_table_arath.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/6/2016.
# Purpose: given a table of TAIR (Arabidopsis) genes and their traits -- as defined either by '0/1' Boolean, or by floating-point value -- along with an Arabidopsis gene_association file from the Gene Ontology Consortium, generate a table usable by FUNC for either hypergeometric or Wilcoxon-rank test.

use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my $trait_table = q{};
my $go_assoc    = q{};
my $no_IEA;
my $help;

my $at_gene = q{};
my $go_term = q{};
my $go_code = q{};
my $trait   = q{};

my $data_ref;

GetOptions ( 'trait_table=s' => \$trait_table,
             'go_assoc=s'    => \$go_assoc,
             'no_IEA'        => \$no_IEA,
             'help'          => \$help,        );

if ($help or (! $trait_table) or (! $go_assoc) ) { 
    die "Format: qual2func_table_arath.pl\n",
        "            -t|--trait_table [TAIR gene trait table]\n",
        "            -g|--go_assoc    [A. thaliana GO-association table from Gene Ontology]\n",
        "            -n|--no_IEA      [ignore IEA annotations]\n",
        "            -h|--help        [print this message]\n",
        "            [print FUNC-usable table to <STDOUT>]\n",
        ;
}

# Sample input lines:
#
# TAIR	locus:2161600	RAN3		GO:0000054	TAIR:Communication:501741973	IBA	PANTHER:PTN000632582	P	AT5G55190	AT5G55190|RAN3|ATRAN3|RAN GTPase 3|MCO15.14|MCO15_14	protein	taxon:3702	20160113	GOC		TAIR:locus:2161600
# 
# TAIR	locus:2031895	AT1G13160		GO:0000055	TAIR:AnalysisReference:501756966	IEA	InterPro:IPR027312	P	AT1G13160	AT1G13160|F3F19.18|F3F19_18	protein	taxon:3702	20160113	TAIR		TAIR:locus:2031895
#
# Note that we are making no effort to select for annotations belonging only to a single ontology, 
#     e.g., 'P -- but we could if we needed to, by a trivial expansion of the code.

open my $GO_ASSOC, '<', $go_assoc or die "Can't open GO-association table $go_assoc: $!";
while (my $input = <$GO_ASSOC>) { 
    chomp $input;
    if ( $input =~ / \A TAIR
                     (?: \t [^\t]*){3}
                     \t (GO[:]\d+)
                     \t [^\t]*
                     \t (I[A-Z]{2})
                     (?: \t [^\t]*){2}
                     \t (AT (?: 1|2|3|4|5|C|M) G\d+) 
                     \t /xms        ) {
        $go_term = $1;
        $go_code = $2;
        $at_gene = $3;
        if ( (! $no_IEA) or ($go_code ne 'IEA') ) { 
            $data_ref->{'gene'}->{$at_gene}->{'GO'}->{$go_term} = 1;
        }
    }
    else { 
        if ( $input !~ /\A # /xms ) { 
            die "Can't parse input line from $go_assoc: $input\n";
        }
    }
}
close $GO_ASSOC or die "Can't close filehandle to GO-association table $go_assoc: $!";

open my $TRAIT_TABLE, '<', $trait_table or die "Can't open trait table $trait_table: $!";
while (my $input = <$TRAIT_TABLE>) { 
    chomp $input;

    # Add '(?: \| \S+ )*' so that I can feed in tables with rich gene names (e.g., 'WBGene00000001|Y110A7A.10|aap-1'), 
    #     rather than just tables with bare IDs like 'WBGene00000001'.

    if ( $input =~ /\A ( (AT (?: 1|2|3|4|5|C|M) G\d+) (?: \| \S+ )* )\t (\S+) \s* \z/xms ) { 
        my $full_gene = $1;
        $at_gene      = $2;
        $trait        = $3;
        if (! looks_like_number($trait) ) { 
            die "Trait $trait does not look like a number; should be boolean 0/1 or floating-point!\n";
        }
        if ( exists $data_ref->{'gene'}->{$at_gene}->{'trait'} ) { 
            die "Gene $at_gene appears to have more than one ascribed trait!\n";
        }
        $data_ref->{'gene'}->{$at_gene}->{'trait'}     = $trait;
        $data_ref->{'gene'}->{$at_gene}->{'full_name'} = $full_gene ;
    }
}
close $TRAIT_TABLE or die "Can't close filehandle to trait table $trait_table: $!";

foreach my $at_gene1 (sort keys %{ $data_ref->{'gene'} } ) {
    # We only care about genes that have *both* traits *and* GO terms; ignore all the others.
    if (      ( exists $data_ref->{'gene'}->{$at_gene1}->{'GO'}    ) 
          and ( exists $data_ref->{'gene'}->{$at_gene1}->{'trait'} ) )  { 

        # The 'full name' can equal the plain WBGene ID -- but it has to exist!
        if (! exists $data_ref->{'gene'}->{$at_gene1}->{'full_name'} ) {
            die "For some reason, can't define full name of WBGene $at_gene1\n";
        }   

        # If we get past all of that, then print the results.
        my $full_name1 = $data_ref->{'gene'}->{$at_gene1}->{'full_name'};
        my $trait1     = $data_ref->{'gene'}->{$at_gene1}->{'trait'};
        foreach my $go_term1 ( sort keys %{ $data_ref->{'gene'}->{$at_gene1}->{'GO'} } ) { 
            print "$full_name1\t$go_term1\t$trait1\n";
        }
    }
}

