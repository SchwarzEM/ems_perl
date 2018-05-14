#!/usr/bin/env perl

# clean_aug_headers.pl -- Erich Schwarz <emsch@its.caltech.edu>, 10/19/2009.
# Purpose: reformat bad getAnnoFasta.pl outputs; die loudly with nonunique name or unparseable header.

use strict;
use warnings;

my %seen = ();

while (my $input = <>) { 
    chomp $input;
    if ( $input !~ /\A > /xms ) { 
        print "$input\n";
    }
    elsif ( $input =~ /\A > ( [^\s;]+  ; ( [^\s;]+ ) .* ) \z /xms ) { 
        my $full_header_text = $1;
        my $seq_id           = $2;
        if ($seen{$seq_id}) { 
            die "Sequence name $seq_id nonredundant!\n";
        }
        print ">$seq_id  $full_header_text\n";
        $seen{$seq_id} = 1;
    }
    # Added this loop as a backstop against really inconsistent exported FASTA headers.
    elsif ( $input =~ /\A > ( [^\s;\.]+  (?:;|\.) ( g [^\s;]+ ) .* ) \z /xms ) {
        my $full_header_text = $1;
        my $seq_id           = $2;
        if ($seen{$seq_id}) {
            die "Sequence name $seq_id nonredundant!\n";
        }
        print ">$seq_id  $full_header_text\n";
        $seen{$seq_id} = 1;
    }
    elsif ( $input =~ /\A > /xms ) { 
        die "Can't parse header: $input\n";
    }
}

