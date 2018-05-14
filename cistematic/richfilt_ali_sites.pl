#!/usr/bin/env perl

# richfilt_ali_sites.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/2/2008.
# Purpose: extract only Ali sites which match enriched motifs.

use strict;
use warnings;

use Getopt::Long;

my $motlist = q{};
my %wanted = ();

GetOptions ( "motifs=s" => \$motlist );

if (! $motlist ) { 
    die "Format: ./richfilt_ali_sites.pl  [-m /--motifs=]MOTLIST  [input site file]\n";
}

open my $MOTLIST, '<', $motlist
    or die "Can't open list of enriched motifs $motlist: $!";
while ( my $input = <$MOTLIST> ) { 
    chomp $input;
    if ($input =~ / \A \s* (\S+) \s* \z /xms ) { 
        my $mot = $1;
        $wanted{$mot} = 1;
    }
}
close $MOTLIST or die "Can't close filehandle to $motlist: $!";

while (my $input = <>) { 
    my $print = 0;
    chomp $input;
    if ( $input !~ / \A \S+ \t \d+ \t \d+ \t \S+ \t [ACGT]+ \z /xms ) {
        warn "Can't parse input line: $input\n";
    }
    if ( $input =~ / \A \S+ \t \d+ \t \d+ \t (\S+) \t [ACGT]+ \z /xms ) {
        my $sourceline = $1;
        my @sources = split /\|/, $sourceline;
        foreach my $source (@sources) { 
            if ( $wanted{$source} ) {
                $print = 1; 
            }
        }
    }
    if ($print) { 
        print "$input\n";
    }
}

