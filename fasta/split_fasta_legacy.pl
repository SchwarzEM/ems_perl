#!/usr/bin/env perl

# split_fasta_legacy.pl -- Erich Schwarz <emsch@its.caltech.edu>, 9/29/2011 -- serious upgrade from 3/13/2008.
# Purpose: split 1+ N-seq. FASTAs (e.g. "supercontigs.fa") into 1+ dirs. of 1-seq. FASTAs; optional suffixes (e.g., '.fa').

use strict;
use warnings;
use Getopt::Long;

my $FASTA;

my $fasta    = q{};
my %seqfiles = ();
my $suffix   = q{};  # Default: no suffix.
my @inputs   = ();
my $date     = q{};

my $directory;
my $help;

GetOptions ( 'inputs=s{,}' => \@inputs,
             'suffix=s'    => \$suffix,
             'dir'         => \$directory,
             'help'        => \$help, );

if ($help or (! @inputs ) ) { 
    die "Format: split_fasta_legacy.pl",
         "--inputs|-i [input files]",
         " --suffixes|-s [suffixes]",
         " --directory|-d [optionally, make an output directory]",
         " --help|-h\n",
         ;
}

if ($directory) { 
    $date = join('.', &get_local_date());
    mkdir $date or die "Can't make directory $date: $!";
    foreach my $infile (@inputs) { 
        mkdir "$date/$infile" or die "Can't make subdirectory $date/$infile: $!";
    }
}

foreach my $infile1 (@inputs) { 
    open my $INFILE, '<', $infile1 or die "Can't open input file: $infile1: $!";
    while (my $input = <$INFILE>) { 
        if ($input =~ /\A > (\S+) /xms) { 
            $fasta = $1;
            if (defined $FASTA) { 
                close $FASTA;
            }
            if ($seqfiles{$fasta}) { 
                warn "Redundant sequence name: $fasta from $infile1\n";
            }
            my $target = $fasta . $suffix;
            if ($directory) { 
                $target = $date. q{/} . $infile1 . q{/} . $target;
            }
            open $FASTA, '>', $target or die "Can't open sequence file $target: $!";
            $seqfiles{$fasta} = 1;
            print $FASTA $input;
        }
        elsif ($fasta) { 
            print $FASTA $input;
        } 
        else { 
            die "Input line does not have a specified FASTA file as destination.\n";
        }
    }
}

sub get_local_date { 
    my @ltime = localtime;
    my @ldate = ( (sprintf ("%04u", ($ltime[5] + 1900)) ),     # $year
                  (sprintf ("%02u", ($ltime[4] + 1))    ),     # $mon
                  (sprintf ("%02u", ($ltime[3] + 0))    ),     # $mday
                  (sprintf ("%02u", ($ltime[2] + 0))    ),     # $hour
                  (sprintf ("%02u", ($ltime[1] + 0))    ),     # $min
                  (sprintf ("%02u", ($ltime[0] + 0))    ), );  # $sec
    return @ldate;
}

