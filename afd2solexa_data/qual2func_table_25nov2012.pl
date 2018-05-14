#!/usr/bin/env perl

# trait2func_table.pl -- Erich Schwarz <emsch@caltech.edu>, 11/25/2012.
# Purpose: given a table of genes and their traits -- as defined either by '0/1' Boolean, or by floating-point value -- along with a gene_association file, generate a table usable by FUNC for either hypergeometric or Wilcoxon-rank test.

use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my $trait_table = q{};
my $go_assoc    = q{};
my $help;

my $gene    = q{};
my $go_term = q{};
my $trait   = q{};

my $data_ref;

GetOptions ( 'trait_table=s' => \$trait_table,
             'go_assoc=s'    => \$go_assoc,
             'help'          => \$help,        );

if ($help or (! $trait_table) or (! $go_assoc) ) { 
    die "Format: qual2func_table.pl\n",
        "            -t|--trait_table [Gene trait table]\n",
        "            -g|--go_assoc    [Gene GO-association table]\n",
        "            -h|--help        [print this message]\n",
        "            [print FUNC-usable table to <STDOUT>]\n",
        ;
}

# Sample input lines:
#
# WBGene00000001   GO:0008340

open my $GO_ASSOC, '<', $go_assoc or die "Can't open GO-association table $go_assoc: $!";
while (my $input = <$GO_ASSOC>) { 
    chomp $input;
    if ( $input =~ / \A (\S+) \t (GO[:]\d+) \z /xms ) { 
        $gene    = $1;
        $go_term = $2;
        $data_ref->{'gene'}->{$gene}->{'GO'}->{$go_term} = 1;
    }
    else { 
        die "From GO-association table $go_assoc, can't parse: $input\n";
    }
}
close $GO_ASSOC or die "Can't close filehandle to GO-association table $go_assoc: $!";

open my $TRAIT_TABLE, '<', $trait_table or die "Can't open trait table $trait_table: $!";
while (my $input = <$TRAIT_TABLE>) { 
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \s* \z/xms ) { 
        $gene    = $1;
        $trait   = $2;
        # We only care about genes that have *both* traits *and* GO terms; no point in recording traits for GO-term-less genes.
        if ( ( exists $data_ref->{'gene'}->{$gene}->{'GO'} ) and (looks_like_number($trait) ) ) { 
            if ( exists $data_ref->{'gene'}->{$gene}->{'trait'} ) { 
                die "Gene $gene has >=2 traits: $data_ref->{'gene'}->{$gene}->{'trait'} and $trait\n";
            }
            $data_ref->{'gene'}->{$gene}->{'trait'} = $trait;
        }
    }
    else { 
        die "From trait table $trait_table, can't parse: $input\n";
    }
}
close $TRAIT_TABLE or die "Can't close filehandle to trait table $trait_table: $!";

foreach my $gene1 (sort keys %{ $data_ref->{'gene'} } ) {
    # Again, we only care about care about genes that have *both* traits *and* GO terms; no point in printing genes without both...
    if ( ( exists $data_ref->{'gene'}->{$gene1}->{'trait'} ) and ( exists $data_ref->{'gene'}->{$gene1}->{'GO'} ) ) { 
        my $trait1 = $data_ref->{'gene'}->{$gene1}->{'trait'};
        foreach my $go_term1 ( sort keys %{ $data_ref->{'gene'}->{$gene1}->{'GO'} } ) { 
            print "$gene1\t$go_term1\t$trait1\n";
        }
    }
}

