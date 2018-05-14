#!/usr/bin/perl

# gff2simplemask.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/18/2008.
# Purpose: softmask simple reps. in FASTA (maybe pre-masked) and GFF2/3 files.

# Caveats:
#
# 1. This script presupposes correct FASTA format, and does nothing to 
#    letter-casing (to allow pre-masking).
# 
# 2. As written on 7/3/2006, this script will fail on any assembly 
#    that marks gaps with 'N' -- i.e., almost all of them!
#    It needs to be revised to take that into account, if it's going 
#    to be used for anything but the totally, completely gap-free C. 
#    elegans chromosomal sequences.

use strict;
use warnings;

my $i              = q{};
my $fasta_file     = q{};
my $gff_file       = q{};
my $input_line     = q{};
my $seq_name       = q{};
my %unmasked_seqs  = ();
my $start_nt       = q{};
my $stop_nt        = q{};
my %mask_ranges    = ();
my $output_line    = q{};
my @unmasked_chars = ();
my @output_lines   = ();

unless ($#ARGV == 1) { 
die 'Format: gff2simplemask.pl ',
    '[FASTA, maybe pre-masked] ',
    '[GFF2/3 file with simple repeats]',
    "\n",
    ; 
}

($fasta_file, $gff_file) = @ARGV;

open my $FASTA, "<", "$fasta_file"
    or die "Can't open FASTA file $fasta_file: $!";

while (<$FASTA>) {
    chomp ($input_line = $_);
    if ($input_line =~ /^>(\S+)/) {
        $seq_name = $1;
        $unmasked_seqs{$seq_name} = q{};
    }
    elsif ($input_line =~ /[a-zA-Z]/) {
        $unmasked_seqs{$seq_name} .= $input_line;
    }
}
close $FASTA;

open my $GFF, "<", "$gff_file"
    or die "Can't open GFF file $gff_file: $!";

my $pattern1 = '(\S+) \s tandem   \s tandem_repeat   \s (\d+) \s (\d+) \s';
my $pattern2 = '(\S+) \s inverted \s inverted_repeat \s (\d+) \s (\d+) \s';

while (<$GFF>) {
    chomp ($input_line = $_);

# sample acceptable input (tab-delimited) lines:
# CHROMOSOME_I    tandem  tandem_repeat   496605  496635  95      .       .       Note "4 copies of 9mer"
# CHROMOSOME_I    inverted        inverted_repeat 7493375 7493484 71      .       .       Note "loop 298, 4 gaps"

    # Simple '##' in regex totally invalidated this conditional!  Use '\#\#'!
    if  ( ($input_line !~ / \A \#\# /xms) and 
          ( ($input_line =~ /\A $pattern1 /xms) or
            ($input_line =~ /\A $pattern2 /xms) ) ) {
        $seq_name = $1;

        # Convert from 1+ DNA counting to 0+ Perl-array counting:
        $start_nt = ($2 - 1);  

        # Ditto:
        $stop_nt  = ($3 - 1);

        unless ($unmasked_seqs{$seq_name} =~ /[A-Z]+/xms) {
            die "Sequence $seq_name with simple repeats ",
                'has no matching FASTA file!',
                "\n",
                ;
        }

        foreach my $nt ($start_nt .. $stop_nt) { 
            $mask_ranges{$seq_name}->{$nt} = 1;
        }
    }
}
close $GFF;

foreach $seq_name (sort keys %unmasked_seqs) { 
    # Rezero this in each loop!
    $output_line = q{}; 

    @unmasked_chars = split(//, $unmasked_seqs{$seq_name});

    $i = 0;
    while ($i < length($unmasked_seqs{$seq_name})) { 
        if ( ( defined $mask_ranges{$seq_name}->{$i}  ) 
             and ( $mask_ranges{$seq_name}->{$i} == 1 ) ) { 
            $unmasked_chars[$i] =~ tr/A-Z/a-z/;
            $output_line .= $unmasked_chars[$i];
        }
        elsif ($unmasked_chars[$i] eq "N") { 
            die 'Hardmasked residue in ',
                'putative un-hardmasked sequence!',
                "\n",
                ; 
        }
        else { 
            $output_line .= $unmasked_chars[$i]; 
        }
        ++$i;
    }

    print ">$seq_name\n";

    @output_lines 
        = unpack("a60" x (length($output_line)/60 + 1), $output_line);

    foreach $output_line (@output_lines) {
          print "$output_line\n";
    }
}

