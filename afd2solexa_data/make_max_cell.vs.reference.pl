#!/usr/bin/env perl

# make_max_cell.vs.reference.pl -- Erich Schwarz <ems394@cornell.edu>, 11/16/2013.
# Purpose: given a set of RPKMs (or FPKMs) of various experimental types, and *one* reference value, compute the highest exp./ref. ratio for each gene.  Typical experimental RPKMs would be from single-cell RT-PCRs of C. elegans using computational pipeline of Schwarz et al. (2012), PubMed ID 22991463; typical reference would be the whole-larval data from that same paper.  Note that the *reference* data (though not the experimental data) needs to have pseudozero values in place of actual 0.00 RPKMs/FPKMs.

use strict;
use warnings;

if (! @ARGV) { 
    die "Format: make_max_cell.vs.reference.pl [list of all experimental data files] [reference data file]\n";
}

# Sample experimental data:
# Gene	CEM_DL_pool
# WBGene00000001|Y110A7A.10|aap-1	0.05
# WBGene00000003|F07C3.7|aat-2	0.05

# Sample reference data:
# Gene	larvae_nz
# WBGene00000001|Y110A7A.10|aap-1	2.35
# WBGene00000003|F07C3.7|aat-2	0.59

my @exp_data = @ARGV;
my $ref_data = pop @exp_data ;

my %gene2maxlc = ();
my %gene2ratio = ();

foreach my $exp_input (@exp_data) { 
    open my $EXP, '<', $exp_input or die "Can't open experimental RPKM input file $exp_input: $!";
    while (my $input = <$EXP>) { 
        chomp $input;
        if ( $input =~ /\A (WBGene\d+\S+) \t (\S+) /xms ) { 
            my $wbgene = $1;
            my $rpkm   = $2;
            if ( ( exists $gene2maxlc{$wbgene} ) and ( $gene2maxlc{$wbgene} < $rpkm ) ) {
                $gene2maxlc{$wbgene} = $rpkm;
            }
            elsif (! exists $gene2maxlc{$wbgene} ) { 
                $gene2maxlc{$wbgene} = $rpkm;
            }
        }
    }
    close $EXP or die "Can't close filehandle to experimental RPKM input file $exp_input: $!";
}

open my $REF, '<', $ref_data or die "Can't open reference input file $ref_data: $!";
while (my $input = <$REF>) { 
    chomp $input;
    if ( $input =~ /\A (WBGene\d+\S+) \t (\S+) /xms ) {
        my $wbgene = $1;
        my $rpkm   = $2;
        if ( (! $gene2maxlc{$wbgene} ) or ( $rpkm == 0 ) ) { 
            die "Can't parse reference data line: $input\n";
        }
        my $ratio = ( $gene2maxlc{$wbgene} / $rpkm);
        $ratio = sprintf("%.2f", $ratio);
        $gene2ratio{$wbgene} = $ratio;
    }
}
close $REF or die "Can't close filehandle to reference input file $ref_data: $!";

my @wbgenes = sort keys %gene2ratio;

print "Gene\t", 'max_experimental/reference', "\n", ;

foreach my $wbgene (@wbgenes) { 
    print "$wbgene\t$gene2ratio{$wbgene}\n";
}

