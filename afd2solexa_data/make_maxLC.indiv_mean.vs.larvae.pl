#!/usr/bin/env perl

use strict;
use warnings;

my @LC_data   = ();
push @LC_data,  $ARGV[0];
push @LC_data,  $ARGV[1];
my $larv_data = $ARGV[2];

my %gene2maxlc = ();
my %gene2ratio = ();

foreach my $lc_input (@LC_data) { 
    open my $LC, '<', $lc_input or die "Can't open LC input file $lc_input: $!";
    while (my $input = <$LC>) { 
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
    close $LC or die "Can't close filehandle to LC input file $lc_input: $!";
}

open my $LARV, '<', $larv_data or die "Can't open larval input file $larv_data: $!";
while (my $input = <$LARV>) { 
    chomp $input;
    if ( $input =~ /\A (WBGene\d+\S+) \t (\S+) /xms ) {
        my $wbgene = $1;
        my $rpkm   = $2;
        if ( (! $gene2maxlc{$wbgene} ) or ( $rpkm == 0 ) ) { 
            die "Can't parse larval data: $input\n";
        }
        my $ratio = ( $gene2maxlc{$wbgene} / $rpkm);
        $ratio = sprintf("%.2f", $ratio);
        $gene2ratio{$wbgene} = $ratio;
    }
}
close $LARV or die "Can't close filehandle to larval input file $larv_data: $!";

my @wbgenes = sort keys %gene2ratio;

print "Gene\t", 'max_wtLC_indiv_mean/larvae_pool', "\n", ;

foreach my $wbgene (@wbgenes) { 
    print "$wbgene\t$gene2ratio{$wbgene}\n";
}


