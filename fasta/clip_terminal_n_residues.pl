#!/usr/bin/env perl

# clip_terminal_n_residues.pl -- Erich Schwarz <ems394@cornell.edu>, 4/26/2015.
# Purpose: given a genome assembly with terminal N or n residues (5' *or* 3'), clip them off.

use strict;
use warnings;
use Getopt::Long;

my @infiles         = ();
my $scaffold        = q{};
my @input_scaffolds = ();

my $data_ref;
my $help;

GetOptions ( 'infiles=s{,}' => \@infiles,
             'help'         => \$help,   );

if ( $help or (! @infiles) ) { 
    die "Format: clip_terminal_n_residues.pl\n",
        "    --infile|-i     <input stream/files>\n",
        "    --help|-h       [print this message]\n",
        ;
}

foreach my $infile (@infiles) { 
    my $INPUT_FILE;
    if ($infile eq '-') {
        # Special case: get the stdin handle
        $INPUT_FILE = *STDIN{IO};
    }
    else {
        # Standard case: open the file
        open $INPUT_FILE, '<', $infile or die "Can't open input file $infile. $!\n";
    }
    while (my $input = <$INPUT_FILE>) { 
        chomp $input;
        if ( $input =~ /\A > (\S+) .*\z/xms ) { 
            $scaffold = $1;
            if ( exists $data_ref->{'scaffold'}->{$scaffold} ) { 
                die "Redundant sequence name: $scaffold\n";
            }
            push @input_scaffolds, $scaffold;

            # Note that $input includes the starting '>' for a FASTA record line.
            $data_ref->{'scaffold'}->{$scaffold}->{'header'} = $input;
        }
        elsif ( $input =~ /\A > /xms ) {
            die "Can't parse input line: $input\n";
        }
        else {
            $data_ref->{'scaffold'}->{$scaffold}->{'seq'} .= $input;
        }
    }
    close $INPUT_FILE or die "Can't close filehandle to input file $infile. $!\n";
}

LOOP: foreach my $scaf (@input_scaffolds) {
    my $seq    = $data_ref->{'scaffold'}->{$scaf}->{'seq'};
    my $header = $data_ref->{'scaffold'}->{$scaf}->{'header'};

    $seq =~ s/\A[Nn]+//;
    $seq =~ s/[Nn]+\z//;
   
    print "$header\n";

    my @output_lines = unpack("a60" x (length($seq)/60 + 1), $seq);
    foreach my $output_line (@output_lines) { 
        if ($output_line =~ /\S/) { 
            print "$output_line\n";
        }
    }
}

