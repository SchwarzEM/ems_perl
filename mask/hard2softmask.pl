#!/usr/bin/perl

# hard2softmask.pl -- Erich Schwarz <emsch@its.caltech.edu>, 7/1/2006; slight modification on 11/28/2013.
# Purpose: make a softmasked FASTA from an unmasked + hardmasked FASTA pair; keep the headers from the unmasked FASTA.

use strict;
use warnings;

my $i = 0;
my $j = 1;

my $unmasked_file   = q{};

my $hardmasked_file = q{};
my $input_line      = q{};

my $seq_name        = q{};
my $header          = q{};

my %seq2header      = ();
my %unmasked_seqs   = ();
my %masked_seqs     = ();

my @unmasked_chars  = ();
my @masked_chars    = ();

my $output_line     = q{};
my @output_lines    = ();

unless ($#ARGV == 1) { die "Format: hard2softmask.pl [unmasked FASTA] [hardmasked FASTA]\n"; }
($unmasked_file, $hardmasked_file) = @ARGV;

open (UNMASKED, "$unmasked_file") or die "Can't open $unmasked_file: $!";
while (<UNMASKED>) { 
    chomp ($input_line = $_);
    if ($input_line =~ /\A > ((\S+).*) \z/xms ) { 
        $header                   = $1;
        $seq_name                 = $2; 
        $unmasked_seqs{$seq_name} = q{};
        $seq2header{$seq_name}    = $header;
    }
    elsif ($input_line =~ /[a-zA-Z]/) { 
        $input_line =~ s/[^a-zA-Z]//g;
        $input_line =~ tr/[a-z]/[A-Z]/;
        $unmasked_seqs{$seq_name} .= $input_line;
    }
}
close UNMASKED;

open (MASKED, "$hardmasked_file") or die "Can't open $hardmasked_file: $!";
while (<MASKED>) {
    chomp ($input_line = $_);
    if ($input_line =~ /\A > (\S+) /xms ) {
        $seq_name = $1;
        unless ($unmasked_seqs{$seq_name} =~ /[A-Z]+/) { 
            die "Masked $seq_name has no unmasked equivalent!\n";
        }
        $masked_seqs{$seq_name} = q{};
    }
    elsif ($input_line =~ /[a-zA-Z]/) {
        $input_line =~ s/[^a-zA-Z]//g;
        $input_line =~ tr/[a-z]/[A-Z]/;
        $masked_seqs{$seq_name} .= $input_line;
    }
}
close MASKED;

foreach my $seq_name1 (sort keys %unmasked_seqs) { 
    $output_line = q{};   # rezero this in each loop!
    unless ( length($masked_seqs{$seq_name1}) == length($unmasked_seqs{$seq_name1}) ) { 
        die "Unmasked $seq_name1 lacks a masked sequence of equal length!\n";
    }
    @unmasked_chars = split(//, $unmasked_seqs{$seq_name1});
    @masked_chars   = split(//, $masked_seqs{$seq_name1});
    $i = 0;
    while ($i < length($unmasked_seqs{$seq_name1})) { 
        if ($unmasked_chars[$i] eq $masked_chars[$i]) { $output_line .= $unmasked_chars[$i]; }

        # This equality test works on any residue, including 'N'.  See consequence below.

        elsif ($masked_chars[$i] eq "N") { 
            $unmasked_chars[$i]  =~ tr/A-Z/a-z/;
            $output_line .= $unmasked_chars[$i];
        }

        elsif ($unmasked_chars[$i] eq "N") { die "Masked residue in putative unmasked sequence!\n"; }

        # This filter works, but there's a really amusing counter-exception. 
        # In supercontig files, it appears that both masked and *un*masked files use 'N'
        #     indiscriminately both to denote masking and to denote gaps between contigs.
        # I found this out running this script on remanei supercontigs (7/17/2006) and finding 
        #     that this script *correctly* parsed that distinction and did *not* die!  
        # Because ... it never ran this test.
        # Instead, it noticed that the 'N' in the unmasked matched an 'N' in the masked,
        #     and silently Did What I Meant!

        elsif ($unmasked_chars[$i] ne $masked_chars[$i]) { 
            $j = $i + 1;
            die "Mismatch of character $j in $seq_name1: $unmasked_chars[$i] versus $masked_chars[$i]!\n";
        }

        ++$i;
    }
    $header = $seq2header{$seq_name1};
    print ">$header\n";
    @output_lines = unpack("a60" x (length($output_line)/60 + 1), $output_line);
    foreach $output_line (@output_lines) {
          print "$output_line\n";
    }
}

