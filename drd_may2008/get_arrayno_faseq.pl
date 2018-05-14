#!/usr/bin/env perl

# get_arrayno_faseq.pl -- Erich Schwarz <emsch@its.caltech.edu>, 5/11/2008.
# Purpose: get Perl-style number of a given seq. in a multiseq. FASTA, for Thornton Scorecons at...
# http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/valdar/scorecons_server.pl

use strict;
use warnings;

my $query_id = $ARGV[0];
my $fasta    = $ARGV[1];

if ( ( $query_id !~ /\w/xms ) or ( $#ARGV != 1) ) { 
    die "Format: ./get_arrayno_faseq.pl [query ID] [FASTA file]\n";
}

my $i = -1;

open my $FASTA, '<', $fasta 
    or die "Can't open FASTA file $fasta: $!";

LOOP:
while (my $input = <$FASTA>) { 
    chomp $input;
    if ( $input =~ /\A > (\S+) /xms ) { 
        my $name = $1;
        $i++;
        if ( $name eq $query_id ) { 
            print "$name is entry $i in FASTA $fasta (first entry being 0).\n";
            last LOOP;
        }
    }
}
 
close $FASTA 
    or die "Failed to close filehandle to FASTA file $fasta: $!";

