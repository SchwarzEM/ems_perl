#!/usr/bin/perl

# uniform_ribo_fasta.pl -- Erich Schwarz <emsch@its.caltech.edu>, 9/27/2007.
# Purpose: change FASTA file (>=1 seqs.) to alph.-ord. headers and ALL-CAPS 60-char/line; make 'T' into 'U'.

use strict;
use warnings;

my $header       = q{};
my $seq_name     = q{};
my %sequences    = ();
my %seq2header   = ();

if ($#ARGV > 0) { 
    die "Error: merging multiple files into one output!\n"; 
} 

while (my $input_line = <>) { 
    chomp $input_line;
    if ($input_line =~ / \A > ( (\S+) .* ) \z /xms) { 
        $header   = $1;
        $seq_name = $2; 
        $sequences{$seq_name}  = q{};
        $seq2header{$seq_name} = $header;
    }
    elsif ($input_line =~ /[a-zA-Z]/) { 
        $input_line =~ s/[^a-zA-Z]//g;
        $input_line =~ tr/[a-z]/[A-Z]/;
        # For some reason I don't quite understand, it's absolutely
        #   critical that this be an 's/T/U/g' and not a 'tr/T/U/':
        $input_line =~ s/T/U/g;
        if ( $input_line =~ /[^ACGU]/ ) { 
            die "Non-RNA residues in: $input_line\n";
        }
        $sequences{$seq_name} .= $input_line;
    }
}

foreach my $seq_name2 (sort keys %sequences) { 
    print ">$seq2header{$seq_name2}\n";
    my @output_lines 
        = grep { $_ =~ /\S/xms } 
          unpack("a60" x (length($sequences{$seq_name2})/60 + 1), 
                 $sequences{$seq_name2});
    foreach my $output_line (@output_lines) { 
          print "$output_line\n";
    }
}

