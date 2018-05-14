#!/usr/bin/env perl

# tbaize_headers.pl -- Erich Schwarz <emsch@its.caltech.edu>, 8/9/2009.
# Purpose: given a FASTA file, convert its header names to the elaborate TBA/multiz-preferred format.

use strict;
use warnings;

foreach my $infile (@ARGV) { 
    my %name2seq = ();
    my $name = q{};
    open my $IN, '<', $infile 
        or die "Can't open input file $infile\n";
    while (my $input = <$IN>) { 
        chomp $input;
        if ( ( $input !~ /\A > /xms ) 
           and ( $input =~ / \S /xms ) ) { 
             if ($name) { 
                 $input =~ s/\s//g;
                 $name2seq{$name} .= $input;
             }
             if (! $name) {
                 die "No sequence name for input line: $input\n";
             }
        }
        if ( $input =~ / \A > (\S+) /xms ) { 
            $name = $1;
        }
        if ( $input =~ / \A > \s* \z /xms ) { 
            die "Aberrant input: $input\n";
        }
    }
    close $IN or die "Can't close filehandle to $infile\n";
    my $outfile = $infile . '.out';
    open my $OUT, '>', $outfile 
        or die "Can't open output file $outfile\n";
    foreach my $outseq (sort keys %name2seq) {
        my $seqlen = length($name2seq{$outseq});
        print $OUT ">$infile:$outseq:1:+:$seqlen\n";
        my @output_lines 
        = unpack( "a60" 
                  x (length($name2seq{$outseq})/60 + 1), 
                  $name2seq{$outseq});
        foreach my $output_line (@output_lines) { 
            if ($output_line =~ /\S/) { 
                print $OUT "$output_line\n";
            }
        }
    }
    close $OUT or die "Can't close filehandle to $outfile\n";
}

