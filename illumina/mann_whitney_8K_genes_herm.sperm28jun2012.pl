#!/usr/bin/env perl

use strict;
use warnings;
use Statistics::Test::WilcoxonRankSum;

my @reinke_genes     = ();
my @non_reinke_genes = ();

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A WBGene\d+\S+ \t (\S+) \t (\S+) (?:\t [^\t]*){44} \t (\S*) \t /xms ) { 
        my $rpkm_l3 = $1;
        my $rpkm_l4 = $2;
        my $reinke  = $3;
        if ( $rpkm_l3 == 0 ) { 
            $rpkm_l3 = 0.01;
        }
        if ( $rpkm_l4 == 0 ) {
            $rpkm_l4 = 0.01;
        }
        my $ratio_l4_l3 = ( $rpkm_l4 / $rpkm_l3 );
        if ( $reinke =~ /herm_sperm/xms ) { 
            push @reinke_genes, $ratio_l4_l3;
        }
        else { 
            push @non_reinke_genes, $ratio_l4_l3;
        }
    }
    else { 
       die "Can't parse input: $input\n";
    }
}

my $wilcox_test = Statistics::Test::WilcoxonRankSum->new();

$wilcox_test->load_data(\@non_reinke_genes, \@reinke_genes);

print "Dataset 1: non-Reinke genes\n";
print "Dataset 2: Reinke genes\n";

my $prob = $wilcox_test->probability();
my $pf    = sprintf '%f', $prob;
print "$pf\n";

my $pstatus = $wilcox_test->probability_status();
print "$pstatus\n";

my $psum = $wilcox_test->summary();
print "$psum\n";

