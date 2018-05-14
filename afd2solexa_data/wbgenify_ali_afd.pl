#!/usr/bin/perl

# wbgenify_ali_afd.pl -- Erich Schwarz <emsch@its.caltech.edu>, 12/24/2007.
# Purpose: map WS170 + WS180 protein-coding gene IDs from headers of wormpep170/180 to Ali AFD data.

use strict;
use warnings;

my %cds2wbgene = ();
my %cds2locus  = ();

unless ($#ARGV == 2) {
    die "Format: ./wbgenify_ali_afd.pl  wormpep170  wormpep180  ali_afd.txt\n";
}

my $input = q{};
my $infile1 = $ARGV[0];
my $infile2 = $ARGV[1];
my $infile3 = $ARGV[2];

open (my $WORMPEP170, "<", "$infile1") 
    or die "Can't open $infile1: $!";
while ($input = <$WORMPEP170>) { 
    chomp $input;
    if ( $input =~ /\A > (\S+) \s+ (.+) \z /xms ) {
        my $cds = $1;
        my $header = $2;
        $cds =~ s/[a-z]+\z//;
        if ($header !~ / (WBGene\d+) /xms ) {
            die "Couldn't parse: $input\n";
        }
        if ($header =~ / (WBGene\d+) /xms ) {
            $cds2wbgene{$cds} = $1;
        }
    }
}
close $WORMPEP170;

open (my $WORMPEP180, "<", "$infile2")
    or die "Can't open $infile2: $!";
while ($input = <$WORMPEP180>) {
    chomp $input;
    if ( $input =~ /\A > (\S+) \s+ (.+) \z /xms ) {
        my $cds = $1;
        my $header = $2;
        $cds =~ s/[a-z]+\z//;
        if ($header !~ / (WBGene\d+) /xms ) {
            die "Couldn't parse: $input\n";
        }
        if ($header =~ / locus:(\S+) /xms ) {
            $cds2locus{$cds} = $1;
        }
    }
}
close $WORMPEP180;


open (my $ALI_FILE, "<", "$infile3") 
    or die "Can't open $infile3: $!";
while ($input = <$ALI_FILE>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \s /xms ) { 
        my $cds = $1;
        if ($cds2wbgene{$cds}) { 
            print "$cds2wbgene{$cds}|";
        }
        if ($cds2locus{$cds}) { 
            print "$cds2locus{$cds}|";
        }
        print "$input\n";
    }
}
close $ALI_FILE;

