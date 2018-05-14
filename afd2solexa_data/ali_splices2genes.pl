#!/usr/bin/perl

# ali_splices2genes.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/19/2008.
# Purpose: given CDS->WBGene list, map Ali-Solexa splices (@ CDS names).

use strict;
use warnings;

my $chromosome = 'I|II|III|IV|V|X';

unless ($#ARGV == 1) { 
    die "Format: ./ali_splices2genes.pl [namefile] [splicefile]\n";
}

my $wormpep    = $ARGV[0];
my $splicefile = $ARGV[1];

my %cds2wbgene = ();
my %wbgene_count = ();

open my $WPEP, "<", "$wormpep" 
    or die "Can't open wormpep file $wormpep\n";

while ( my $input = <$WPEP> ) { 
    chomp $input;

    # Sample input lines:
    # >4R79.1b CE39659 WBGene00003525 locus:nas-6 ... 
    # >4R79.2 CE19650 WBGene00007067 ...

    if ( $input =~ / \A >
                     ( \w+\.\w+ )   # $1; CDS
                     \s .+?
                     ( WBGene\d+ )  # $2; WBGene ID no.
                     \s
                   /xms ) {
        my ($cds, $wbgene);
        ($cds, $wbgene) = ($1, $2);
        $cds =~ s/[a-z]\z//;
        $cds2wbgene{$cds} = $wbgene;
    }
}

close $WPEP;

open my $SPLICES, "<", "$splicefile" 
    or die "Can't open splice file $splicefile\n";

while ( my $input = <$SPLICES> ) { 
    chomp $input;

    # Sample input line:
    # chrIII  4753789 4753886 B0393.1-7F1 ...

    if ( $input =~ / chr(?: $chromosome )  # not $1 !
                     \s+ \d+ \s+ \d+ \s+
                     ( \w+\.\w+ )          # e.g. "B0393.1"
#                    \- \S+
                   /xms ) { 
        my $cds;
        $cds = $1;
        $cds =~ s/[a-z]\z//;
        if ( exists $cds2wbgene{$cds} ) { 
            my $wbgene = $cds2wbgene{$cds};
            $wbgene_count{$wbgene}++;
        }
        if (! exists $cds2wbgene{$cds} ) { 
            warn "Couldn't find WBGene name for CDS name $cds\n";
        }
    }
}

close $SPLICES;

my @genes = sort 
    { $wbgene_count{$b} <=> $wbgene_count{$a} } keys %wbgene_count;

foreach my $gene (@genes) { 
    print "$gene\t$wbgene_count{$gene}\n";
}

