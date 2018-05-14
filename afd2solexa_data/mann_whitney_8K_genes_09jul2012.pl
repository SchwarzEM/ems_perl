#!/usr/bin/env perl

use strict;
use warnings;
use Statistics::Test::WilcoxonRankSum;

my @reinke_genes     = ();
my @non_reinke_genes = ();

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A WBGene\d+\S+ \t (\S+) \t (\S+) \t (\S*) /xms ) { 
        my $rpkm_l3     = $1;
        my $rpkm_nonl3 = $2;
        my $reinke      = $3;
        if ( $rpkm_l3 == 0 ) { 
            $rpkm_l3 = 0.01;
        }
        if ( $rpkm_nonl3 == 0 ) {
            $rpkm_nonl3 = 0.01;
        }
        my $ratio_nonl3_l3 = ( $rpkm_nonl3 / $rpkm_l3 );
        if ( $reinke =~ /sperm/xms ) { 
            push @reinke_genes, $ratio_nonl3_l3;
        }
        else { 
            push @non_reinke_genes, $ratio_nonl3_l3;
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

