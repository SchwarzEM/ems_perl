#!/usr/bin/env perl

# uniq_GFF_seqlines.pl -- Erich Schwarz <emsch@its.caltech.edu>, 7/11/2008.
# Purpose: have one sequence line per segment in GFF; act as filter.

use strict;
use warnings;

my $sequence     = q{};
my %stored_lines = ();
my %do_not_store = ();

while (my $input = <>) { 
    chomp $input;

    # Store ## lines w/o immediate printing.
    # E.g.: ##sequence-region Contig1000 1 15788
    if (  ( $input =~ / \A 
                        [#]{2}
                        sequence-region 
                        \s+ 
                        (.+?) 
                        \s+ 
                        1 
                        \s 
                        \d+ /xms ) 
          and (! $do_not_store{$1}) ) { 
        $sequence = $1;
        $stored_lines{$sequence} = $input;
    }

    # Print Sequence GFF lines; unstore ##; block more storage of them.
    # E.g.: Contig1000  .  Sequence  1  15788  .  +  .  Sequence "Contig1000"
    elsif ( $input =~ / \A 
                        ([^\t]+) 
                        \t 
                        \. 
                        \t 
                        Sequence 
                        \t 
                        1 
                        \t 
                        \d+ /xms ) { 
        $sequence = $1;
        if ($stored_lines{$sequence}) { 
            delete $stored_lines{$sequence};
            $do_not_store{$sequence} = 1;
        }
        print "$input\n";
    }

    # Just print everything else.
    else { 
        print "$input\n";
    }
}

foreach my $unlined (sort keys %stored_lines) { 
    print "$stored_lines{$unlined}\n";
}

