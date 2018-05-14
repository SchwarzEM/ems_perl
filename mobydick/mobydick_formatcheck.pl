#!/usr/bin/env perl

# mobydick_formatcheck.pl -- Erich Schwarz <emsch@its.caltech.edu>, 10/21/2008.
# Purpose: validate what I think is mobydick sequence file format.
 
use strict;
use warnings;

my $lines_scanned  = 0;
my $lines_accepted = 0;
my $comment_lines  = 0;

while (my $input = <>) { 
    chomp $input;
    $lines_scanned++;
    if ( ( $input !~ /\A \# /xms ) 
         and ($input !~ /\A \S+ : \d+ : \- \d+ : [acgt]+ : \z /xms ) ) { 
        warn "Non-mobydick format line: $input\n";
    }
    if ( $input =~ /\A \# /xms ) { 
        $lines_accepted++;
        $comment_lines++;
    }
    if ($input =~ /\A \S+ : (\d+) : \- (\d+) : ([acgt]+) : \z /xms ) { 
        my $len1 = $1;
        my $len2 = $2;
        my $seq  = $3;
        my $seqlen = length($seq);
        if ( $seqlen != $len1 ) { 
            warn "Non-matching sequence lengths in: $input\n";
            warn "Claimed $len1 versus real $seqlen\n";
        }
        if ( $seqlen == $len1 ) {
            $lines_accepted++;
        }
    }
}

print "SUCCESS: $lines_accepted out of $lines_scanned lines";
if ($comment_lines) { 
    print " (with $comment_lines comment lines)";
}
print " were scanned without format errors.\n";

