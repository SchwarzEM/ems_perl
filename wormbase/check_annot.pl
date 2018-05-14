#!/usr/bin/env perl

# check_annot.pl -- Erich Schwarz <emsch@its.caltech.edu>, 8/8/08
# Purpose: scan the first five words of concise description for a gene.

use strict;
use warnings;

my $EXCERPT_LEN = 5;
my %keywords   = ();
my @to_read    = ();

foreach my $word (@ARGV) { 
    if (! -e $word) { 
        $keywords{$word} = 1;
    }
    if (-e $word) { 
        push @to_read, $word;
    }
}

foreach my $infile (@to_read) { 
    open my $FILE, '<', $infile or die "Can't open file $infile: $!";
    while (my $input = <$FILE>) { 
        chomp $input;
        if ( $input =~ / \A Concise_description \s+ \"
                         ( ( (?: \S+ \s+ ){$EXCERPT_LEN} ) .* ) \" /xms ) { 
            my $inline  = $1;
            my $excerpt = $2;

            my @words = split /\s+/, $excerpt;
            foreach my $word (@words) { 
                if ( $keywords{$word} ) { 
                    print "\n$inline\n\n";
                }
            }
        }
    }
    close $FILE or die "Can't close filehandle of $infile: $!";
}

