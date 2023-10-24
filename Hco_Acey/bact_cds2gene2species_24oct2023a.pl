#!/usr/bin/env perl

use strict;
use warnings;
use strict;

my $infile  = q{};
my $species = q{};

$infile  = $ARGV[0] if $ARGV[0];
$species = $ARGV[1] if $ARGV[1];

if ( (! $infile ) or (! $species ) ) {
    die "Format: bact_cds2gene2species_24oct2023a.pl [infile] [species] > [cds2gene2species]\n";
}

open my $INFILE, '<', $infile;
while ( my $input = <$INFILE> ) {
    chomp $input;
    if ( $input =~ /\A [>] (\S+ \| \S+ \| (\S+) \|) \s/xms ) {
        my $cds  = $1;
        my $gene = $2;
        $gene    = $species . '-' . $gene;
        print "$cds\t$gene\t$species\n";
    }
    elsif ( $input =~ /\A [>] ((\S+) \| \S+ \|) \s/xms ) {
        my $cds  = $1;
        my $gene = $2;              
        $gene    = $species . '-' . $gene;
        print "$cds\t$gene\t$species\n";
    }
    elsif ( $input =~ /\A [>] /xms ) {
        die "From infile $infile, cannot parse input header line: $input\n";
    }
}
close $INFILE;

# >PAU_00004|ID:3711956|gyrB| DNA gyrase subunit B [Photorhabdus asymbiotica ATCC43949]
# >PHAv3_0005|ID:3710932| protein of unknown function [Photorhabdus asymbiotica ATCC43949]



sub safename {
    my $filename = $_[0];
    my $orig_filename = $filename;
    if (-e $orig_filename) {
        my $suffix1 = 1;
        $filename = $filename . ".$suffix1";
        while (-e $filename) {
            $suffix1++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix1";
        }
    }
    return $filename;
}

