#!/usr/bin/env perl

# map_rscouts_to_wpep.pl -- Erich Schwarz <emsch@its.caltech.edu>, 8/24/2008.
# Purpose: given BlastX/RepeatScout/wormpep, map any hits to full headers.

use strict;
use warnings;

my $wpep        = $ARGV[0];
my $reps_blastx = $ARGV[1];

my $cds        = q{};
my %cds2header = ();

my $repeat     = q{};
my %repeat2cds = ();
my %repeat2sim = ();

my $scanning   = 0;

open my $WPEP, '<', $wpep 
    or die "Can't open wormpep (header) file: $!";
while ( my $input = <$WPEP> ) { 
    chomp $input;
    if ( $input =~ /\A > ( (\S+) .*) \z/xms ) { 
        $cds    = $2;
        my $header = $1;
        $cds2header{$cds} = $header;
    }
}
close $WPEP or die "Can't close filehandle to $wpep: $!";

open my $REPS_BLASTX, '<', $reps_blastx 
    or die "Can't open repeats.blastx.wormpep file $reps_blastx: $!";
LOOP:
while ( my $input = <$REPS_BLASTX> ) { 
    chomp $input;
    if ( $input =~ / \A \s* Query= \s+ (.+) \s* \z/xms) {
        $repeat = $1;
        $cds      = q{};
        $scanning = 0;
    }
    if ( $input =~ /\* \s No \s hits \s found \s \*/xms ) { 
        $repeat   = q{};
        $cds      = q{};
        $scanning = 0;
        next LOOP;
    }
    if ( $input =~ / Sequences 
                     \s producing 
                     \s significant 
                     \s alignments /xms ) { 
        $scanning = 1;
    }
    if ( ( $scanning == 1 ) and ( $input =~ / \A > (\S+) /xms ) ) { 
        $cds = $1;
        $repeat2cds{$repeat} = $cds;
        $scanning = 2;
    }
    if ( ( $scanning == 2 ) 
         and ( $input =~ / Identities .+ Positives /xms ) ) { 
        $input =~ s/\A\s+|\s+\z//g;
        $repeat2sim{$repeat} = $input;
        $repeat   = q{};
        $cds      = q{};
        $scanning = 0;
        next LOOP;
    }
}
close $REPS_BLASTX or die "Can't close filehandle to $reps_blastx: $!";

foreach my $positive (sort keys %repeat2sim) { 
    my $match_cds    = $repeat2cds{$positive};
    my $match_sim    = $repeat2sim{$positive};
    my $match_header = $cds2header{$match_cds};
    my $output = "$positive\t$match_header\t$match_sim";
    $output =~ s/[ ]{2,}/ /g;
    print "$output\n";
}

