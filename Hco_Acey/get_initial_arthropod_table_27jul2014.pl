#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my %pref2sp = (
    Acpis => 'Acyrthosiphon pisum',
    Aeaeg => 'Aedes aegypti',
    Angam => 'Anopheles gambiae',
    Atcep => 'Atta cephalotes',
    Bomor => 'Bombyx mori',
    Cuqui => 'Culex quinquefasciatus',
    Daple => 'Danaus plexippus',
    Dapul => 'Daphnia pulex',
    Depon => 'Dendroctonus ponderosae',
    Hemel => 'Heliconius melpomene',
    Ixsca => 'Ixodes scapularis',
    Mesca => 'Megaselia scalaris',
    Navit => 'Nasonia vitripennis',
    Rhpro => 'Rhodnius prolixus',
    Soinv => 'Solenopsis invicta',
    Stmar => 'Strigamia maritima',
    Teurt => 'Tetranychus urticae',
    Trcas => 'Tribolium castaneum',
);

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A > (([A-Z][a-z]+) [_] (\S+)) \b .* \b gene:(\S+) /xms ) { 
        my $seqname = $1;
        my $prefix  = $2;
        my $accname = $3;
        my $gene    = $4;
        if (! exists $pref2sp{$prefix} ) {
            die "Can't parse species prefix in: $input\n";
        }
        my $species = $pref2sp{$prefix};
        $seqname = 'Arth_' . $seqname;
        print "$seqname\tArthropods\t$species\t$accname\t$gene\n";
    }
    elsif ( $input =~ /\A > /xms ) { 
        die "Can't parse FASTA sequence header: $input\n";
    }
}

