#!/usr/bin/env perl

# summarize_2011_RNAi_hits.pl -- Erich Schwarz <emsch@its.caltech.edu>, 6/5/2011.
# Purpose: get a clean table listing all genes given RNAi, and whether they were 'WT_RNAi' or 'nonwt_RNAi'.

use strict;
use warnings;

my $wbgene      = q{};
my %phenogenes  = ();
my %seen        = ();
my $header_seen = 0;

open my $PHENO_LIST, '<', $ARGV[0] or die "Can't open list of genes with RNAi phenotypes, $ARGV[0]: $!";
while (my $input = <$PHENO_LIST>) { 
    chomp $input;
    if ( $input =~ /\A ( WBGene \d+ \| \S+ ) \s* /xms ) { 
        $wbgene = $1;
        $phenogenes{$wbgene} = 1;
    }
    else { 
        warn "Can't parse input from list of genes with RNAi phenotypes: $input\n";
    }
}
close $PHENO_LIST or die "Can't close filehandle to list of genes with RNAi phenotypes, $ARGV[0]: $!";

print "Gene\tRNAi_screen\n";

open my $RNAI_LIST, '<', $ARGV[1] or die "Can't open list of genes tested with RNAi, $ARGV[1]: $!";

while (my $input = <$RNAI_LIST>) { 
    chomp $input;
    # WBGene00015419|C04E6.2|srsx-18 etc.
    if ( $input =~ /\A ( WBGene \d+ \| \S+ ) \s* /xms ) { 
        $wbgene = $1;
        if (! $seen{$wbgene} ) { 
            print "$wbgene\t";
            if ( exists $phenogenes{$wbgene} ) { 
                print "nonwt_RNAi\n";
            }
            else { 
                print "WT_RNAi\n";
            }
        }
        $seen{$wbgene} = 1;
    }
    else { 
        warn "Can't parse input line from list of genes tested with RNAi: $input\n";
    }
}

close $RNAI_LIST or die "Can't close filehandle to list of genes tested with RNAi, $ARGV[1]: $!";

