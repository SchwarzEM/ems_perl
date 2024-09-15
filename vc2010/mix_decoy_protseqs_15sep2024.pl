#!/usr/bin/env perl

# mix_decoy_protseqs_15sep2024.pl -- Erich Schwarz <ems394@cornell.edu>, 9/15/2024.
# Purpose: interleave a target proteome with a decoy proteome for subsequent mass. spec. analysis of the proteome.

use strict;
use warnings;
use autodie;

use Getopt::Long;

my $target       = q{};
my $decoy        = q{};
my $prefix       = 'decoy_';
my $seq_name     = q{};
my %sequences    = ();
my @target_names = ();
my $help;

GetOptions ( 'target=s' => \$target,
             'decoy=s'  => \$decoy,
             'prefix:s' => \$prefix,
             'help'     => \$help, ); 

if ( $help or (! $target ) or (! $decoy ) ) {
    print "Format: mix_decoy_protseqs_15sep2024.pl\n",
          "            --target|-t  [input target FASTA file]\n",
          "            --decoy|-d   [input decoy FASTA file]\n",
          "            --prefix|-p  [prefix for decoy seqs.; default is \"decoy_\"]\n",
          "            --help|-h    [print this help message\n",
          ;
    exit;
}

open my $DECOY, '<', $decoy;
while (my $input_line = <$DECOY>) {
    chomp $input_line;
    if ($input_line =~ /\A [>] (\S+) /xms) {
        $seq_name = $1;
    }
    elsif ( $input_line =~ /\A > /xms ) {
        die "In decoy FASTA $decoy with last seen sequence $seq_name, aberrant input line: \"$input_line\"\n";
    }
    elsif ( $input_line =~ /[a-zA-Z]/xms ) {
        $sequences{$seq_name} .= $input_line;
    }
}
close $DECOY;

open my $TARGET, '<', $target;
while (my $input_line = <$TARGET>) { 
    chomp $input_line;
    if ($input_line =~ /\A [>] (\S+) /xms) { 
        $seq_name = $1;
        my $decoy_name = "$prefix$seq_name";
        if (! exists $sequences{$decoy_name} ) {
            die "In target FASTA $target, cannot match sequence $seq_name to the decoy sequence $decoy_name\n";
        }
        $sequences{$seq_name} = q{};
        push @target_names, $seq_name;
    }
    elsif ( $input_line =~ /\A > /xms ) { 
        die "In target FASTA $target with last seen sequence $seq_name, aberrant input line: \"$input_line\"\n";
    }
    elsif ( $input_line =~ /[a-zA-Z]/xms ) { 
        $sequences{$seq_name} .= $input_line;
    }
}
close $TARGET;

foreach my $target_name (@target_names) { 
    if (! exists $sequences{$target_name} ) {
        die "Before interleaving target sequence $target_name, cannot find its sequence\n";
    }

    my $decoy_name = "$prefix$target_name"; 
    if (! exists $sequences{$decoy_name} ) {
        die "Before interleaving target sequence $target_name, cannot find sequence for decoy sequence $decoy_name\n";
    }

    # Print a TARGET sequence.
    print ">$target_name\n";
    my @output_target_lines = unpack("a60" x (length($sequences{$target_name})/60 + 1), $sequences{$target_name});
    foreach my $output_target_line (@output_target_lines) { 
        if ($output_target_line =~ /\S/) { 
            print "$output_target_line\n";
        }
    }

    # Then, print its corresponding DECOY sequence.
    print ">$decoy_name\n";
    my @output_decoy_lines = unpack("a60" x (length($sequences{$decoy_name})/60 + 1), $sequences{$decoy_name});
    foreach my $output_decoy_line (@output_decoy_lines) {
        if ($output_decoy_line =~ /\S/) {
            print "$output_decoy_line\n";
        }
    }
}
