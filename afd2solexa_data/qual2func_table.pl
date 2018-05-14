#!/usr/bin/env perl

# trait2func_table.pl -- Erich Schwarz <emsch@its.caltech.edu>, 12/3/2010
# Purpose: given a table of WBGenes and their traits -- as defined either by '0/1' Boolean, or by floating-point value -- along with a WormBase gene_association file, generate a table usable by FUNC for either hypergeometric or Wilcoxon-rank test.

use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my $trait_table = q{};
my $go_assoc    = q{};
my $no_IEA;
my $help;

my $wb_gene = q{};
my $go_term = q{};
my $go_code = q{};
my $trait   = q{};

my $data_ref;

GetOptions ( 'trait_table=s' => \$trait_table,
             'go_assoc=s'    => \$go_assoc,
             'no_IEA'        => \$no_IEA,
             'help'          => \$help,        );

if ($help or (! $trait_table) or (! $go_assoc) ) { 
    die "Format: qual2func_table.pl\n",
        "            -t|--trait_table [WBGene trait table]\n",
        "            -g|--go_assoc    [WB GO-association table]\n",
        "            -n|--no_IEA      [ignore IEA annotations]\n",
        "            -h|--help        [print this message]\n",
        "            [print FUNC-usable table to <STDOUT>]\n",
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

open my $TRAIT_TABLE, '<', $trait_table or die "Can't open trait table $trait_table: $!";
while (my $input = <$TRAIT_TABLE>) { 
    chomp $input;

    # Add '(?: \| \S+ )*' so that I can feed in tables with rich gene names (e.g., 'WBGene00000001|Y110A7A.10|aap-1'), 
    #     rather than just tables with bare IDs like 'WBGene00000001'.

    if ( $input =~ /\A ( (WBGene\d+) \S*? ) \t (\S+) \s* \z/xms ) { 
        my $full_gene = $1;
        $wb_gene      = $2;
        $trait        = $3;
        if (! looks_like_number($trait) ) { 
            die "Trait $trait does not look like a number; should be boolean 0/1 or floating-point!\n";
        }
        if ( exists $data_ref->{'gene'}->{$wb_gene}->{'trait'} ) { 
            die "Gene $wb_gene appears to have more than one ascribed trait!\n";
        }
        $data_ref->{'gene'}->{$wb_gene}->{'trait'}     = $trait;
        $data_ref->{'gene'}->{$wb_gene}->{'full_name'} = $full_gene ;
    }
}
close $TRAIT_TABLE or die "Can't close filehandle to trait table $trait_table: $!";

foreach my $wb_gene1 (sort keys %{ $data_ref->{'gene'} } ) {
    # We only care about genes that have *both* traits *and* GO terms; ignore all the others.
    if (      ( exists $data_ref->{'gene'}->{$wb_gene1}->{'GO'}    ) 
          and ( exists $data_ref->{'gene'}->{$wb_gene1}->{'trait'} ) )  { 

        # The 'full name' can equal the plain WBGene ID -- but it has to exist!
        if (! exists $data_ref->{'gene'}->{$wb_gene1}->{'full_name'} ) {
            die "For some reason, can't define full name of WBGene $wb_gene1\n";
        }   

        # If we get past all of that, then print the results.
        my $full_name1 = $data_ref->{'gene'}->{$wb_gene1}->{'full_name'};
        my $trait1     = $data_ref->{'gene'}->{$wb_gene1}->{'trait'};
        foreach my $go_term1 ( sort keys %{ $data_ref->{'gene'}->{$wb_gene1}->{'GO'} } ) { 
            print "$full_name1\t$go_term1\t$trait1\n";
        }
    }
}

