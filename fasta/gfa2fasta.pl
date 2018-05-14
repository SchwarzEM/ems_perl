#!/usr/bin/env perl

# gfa2fasta.pl -- Erich Schwarz <ems394@cornell.edu>, 12/7/2016.
# Purpose: convert miniasm graph output files (*.gfa) to more generally usable FASTA files.

use strict;
use warnings;
use Getopt::Long;

my $infile       = q{};
my $output_line  = q{};
my @output_lines = ();
my $seq_name     = q{};
my %seqs2headers = ();
my %sequences    = ();
my $allcaps      = 0;
my $fix2residue  = q{};
my @censor       = ();
my $keep_order   = 0;
my @orig_names   = ();
my @sorted_names = ();
my $help;

GetOptions ( 'infile=s'      => \$infile,
             'allcaps'       => \$allcaps,
             'fix_residue=s' => \$fix2residue,
             'censor=s{,}'   => \@censor,
             'keep_order'    => \$keep_order,
             'help'          => \$help, ); 

if ($help or (! $infile) ) {
    print "Format: gfa2fasta.pl\n",
          "            --infile|-i             [input GFA file or '-']\n",
          "            --censor|-c             [silently delete one or more specified residues, e.g, '-']\n",
          "            --fix_residue|-f <ARG>  [correct bad, uncensored residues to a single residue -- usually 'X' or 'N']\n",
          "            --keep_order|-k         [keep original order of sequences in input]\n",
          "            --allcaps|-a\n",
          "            --help|-h               [print this message]\n",
          ;
    exit;
}

# Accept either a stream from '-' or a standard file.
my $INPUT_FILE;
if ($infile eq '-') {
    # Special case: get the stdin handle
    $INPUT_FILE = *STDIN{IO};
}
else {
    # Standard case: open the file
    open $INPUT_FILE, '<', $infile or die "Can't open input file $infile. $!\n";
}

while (my $input_line = <$INPUT_FILE>) { 
    chomp $input_line;
    if ($input_line =~ /\A [S] \t (\S+) \t (\S+) \t (\S+) \z/xms ) { 
        $seq_name = $1;
        my $sequence = $2;
        my $comment  = $3;

        $seqs2headers{$seq_name} = "$seq_name  $comment";
        $sequences{$seq_name} = $sequence;
        if ($keep_order) {
            push @orig_names, $seq_name;
        }
    }
    elsif ( $input_line =~ /\A [S] \t /xms ) { 
        die "Aberrant input line: \"$input_line\"\n";
    }
}
close $INPUT_FILE or die "Can't close filehandle to input file $infile. $!\n";

if ($keep_order) {
    @sorted_names = @orig_names ;
}
else {
    @sorted_names = sort keys %sequences;
}

foreach my $seq_name2 (@sorted_names) { 
    # Do all editing and error-checking with final raw sequence -- this makes catching internal stop codons much easier!
    # Trim terminal stop codons.
    $sequences{$seq_name2} =~ s/\*\z//;

    # Die loudly if given an internal stop codon.
    if ( $sequences{$seq_name2} =~ /\*/xms ) { 
        if (! $fix2residue) { 
            die "FATAL: $seq_name2 contains an internal stop codon!\n";
        }
        if ($fix2residue) { 
            warn "WARNING: $seq_name2 contains an internal stop codon!\n";
        }
    }

    # Remove all whitespace.
    $sequences{$seq_name2} =~ s/\s//g;

    # Optionally, censor any given residues:
    if (@censor) { 
        foreach my $banned (@censor) { 
            $sequences{$seq_name2} =~ s/$banned//g;
        }
    }

    # Optionally, change any weird characters to an acceptable residue (typically 'X' or 'N'):
    if ($fix2residue) { 
        $sequences{$seq_name2} =~ s/[^ACGTNacgtn]/$fix2residue/g;
    }

    # Die loudly if any weird characters remain in the sequence.
    my $badchar = '[unclear]';
    if ( $sequences{$seq_name2} =~ /([^a-zA-Z])/xms ) { 
        $badchar = $1;
        warn "$seq_name2 contains at least one unacceptable character: \"$badchar\"!\n"; 
    }

    # Optionally, write the sequence in GREAT RUNES.
    if ($allcaps) {
        $sequences{$seq_name2} =~ tr/[a-z]/[A-Z]/;
    }

    # At last, print the sequence.
    print ">$seqs2headers{$seq_name2}\n";
    @output_lines 
        = unpack("a60" x (length($sequences{$seq_name2})/60 + 1), $sequences{$seq_name2});
    foreach $output_line (@output_lines) { 
        if ($output_line =~ /\S/) { 
            print "$output_line\n";
        }
    }
}

