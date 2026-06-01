#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $ident_count     = 0;
my $non_ident_count = 0;

while ( my $input = <> ) {
    chomp $input;
    # Sample: ChrII     AUGUSTUS        mRNA    9694    15025   0.19    -       .       ID=CBR03048;Parent=CBR03048
    if ( $input =~ / mRNA /xms ) {
        if ( $input =~ /\A \S+ \t \S+ \t mRNA \t \d+ \t \d+ \t \S+ \t \S+ \t \S+ \t ID=(\S+)[;]Parent=(\S+) \z/xms ) {
            my $tx   = $1;
            my $gene = $2;
            if ( $tx eq $gene ) {
                $ident_count++;
            }
            else {
                $non_ident_count++;
                warn "Instance of non-identical names: gene $gene / transcript $tx\n";
            }
        }
        else {
            warn "Dubious line: $input\n";
        }
    }
}

print "Transcripts with names identical to their genes:    $ident_count\n";
print "Transcripts with names different from their genes:  $non_ident_count\n";

