#!/usr/bin/perl

# weeded_fasta.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/5/2007.
# Purpose: derivative of uniform_fasta.pl that censors zero-length, no M1, or int-stopped proteins.

use strict;
use warnings;
use Getopt::Long;

my $output_line  = q{};
my @output_lines = ();
my $seq_name     = q{};
my %seqs2headers = ();
my %sequences    = ();
my $allcaps      = 0;

GetOptions ("allcaps" => \$allcaps); 

while (my $input_line = <>) { 
    chomp $input_line;
    if ($input_line =~ /\A > ( (\S+) .*) /xms) { 
        $seq_name = $2; 
        $seqs2headers{$seq_name} = $1;
        $sequences{$seq_name} = q{};
    }
    elsif ( $input_line =~ /[a-zA-Z]/xms ) { 
        $input_line =~ s/[^a-zA-Z]//g;
        if ($allcaps) {
            $input_line =~ tr/[a-z]/[A-Z]/;
        }
        $sequences{$seq_name} .= $input_line;
    }
}

foreach my $seq_name2 (sort keys %sequences) { 
    # Require: starting methionine; non-empty sequences; no internal stop codons.
    if ( ( $sequences{$seq_name2} =~ /[a-zA-Z]/xms ) 
       and ( $sequences{$seq_name2} !~ /[a-zA-Z]\*[a-zA-Z]/xms ) 
       and ( $sequences{$seq_name2} =~ / \A (m|M) /xms ) ) { 
        print ">$seqs2headers{$seq_name2}\n";
        @output_lines 
            = unpack("a60" x (length($sequences{$seq_name2})/60 + 1), $sequences{$seq_name2});
        foreach $output_line (@output_lines) { 
            if ($output_line =~ /\S/) { 
                print "$output_line\n";
            }
        }
    }
}

