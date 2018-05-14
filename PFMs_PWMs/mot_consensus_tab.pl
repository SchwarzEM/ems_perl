#!/usr/bin/env perl

# mot_consensus_tab.pl -- Erich Schwarz <emsch@its.caltech.edu>, 5/24/2010.
# Purpose: given a MEME motif text file, make a TSV of motif names and consenses.

use strict;
use warnings;
use Getopt::Long;

my $prefix     = q{};
my $motif      = q{};
my $read       = 0;
my $consensus  = q{};
my @motif_list = ();
my %mot2cons   = ();

my %chars2cons = ( AC   => 'M',
                   AG   => 'R',
                   AT   => 'W',
                   CT   => 'Y',
                   CG   => 'S',
                   GT   => 'K',
                   ACT  => 'H',
                   ACG  => 'V',
                   AGT  => 'D',
                   CGT  => 'B',
                   ACGT => 'N', );

GetOptions ( 'prefix:s' => \$prefix, );

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A \s+ Motif \s+ \b (.+) \b \s+ regular \s+ expression \s* \z /xms ) { 
        $motif     = $1;
        $read      = 0;
        $consensus = q{};
    }
    if ($motif and ( $input =~ /\A [\-]+ \s* \z /xms ) ) { 
        $read = 1;
    }
    if ( $motif and $read and ( $input =~ /\A \s* ([ACGT\[] [ACGT\[\]]*) \s* \z /xms ) ) { 
        $consensus = $1;
        $consensus = regex2cons($consensus);
        if (exists $mot2cons{$motif} ) { 
            die "Redundant consenses listed for motif $motif!\n";
        }
        push @motif_list, $motif;
        $mot2cons{$motif} = $consensus;
        $motif     = q{};
        $read      = 0;
        $consensus = q{};
    }
}

foreach my $mot2 (@motif_list) {
    print "$prefix$mot2\t$mot2cons{$mot2}\n";
}

sub regex2cons { 
    my $_input           = $_[0];
    my $_output          = q{};
    my @_input_chars     = split //, $_input;
    my @_consensus_chars = ();
    my $_getting_cons    = 0;
    foreach my $_char (@_input_chars) { 
        if ( (! $_getting_cons ) and ( $_char =~ /\A [ACGT] \z/xms ) ) { 
            $_output .= $_char;
        }
        if ( (! $_getting_cons ) and ( $_char =~ /\A \[ \z/xms ) ) { 
            $_getting_cons = 1;
        }
        if ( ( $_getting_cons ) and ( $_char =~ /\A [ACGT] \z/xms ) ) {
            push @_consensus_chars, $_char;
        }
        if ( ( $_getting_cons ) and ( $_char =~ /\A \] \z/xms ) ) {
            $_getting_cons = 0;
            @_consensus_chars = sort @_consensus_chars;
            my $_cons_chars = join q{}, @_consensus_chars;
            if (! exists $chars2cons{$_cons_chars} ) { 
                die "Can't generate single consensus symbol for residues \"$_cons_chars\"!\n";
            }
            $_output .= $chars2cons{$_cons_chars};
            @_consensus_chars = ();
        }
    }
    return $_output;
}

