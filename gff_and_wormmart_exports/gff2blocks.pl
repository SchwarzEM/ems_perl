#!/usr/bin/env perl

# gff2blocks.pl -- Erich Schwarz <emsch@caltech.edu>, 10/20/2011.
# Purpose: given a GFF, print simple block lines ("seqname\tstart_nt\tend_nt\n"); options to enforce ascending order of nt sites, or to strip prefixes such as 'chr' from seqnames.

use strict;
use warnings;
use Getopt::Long;

my @input_files = ();
my $ascending;
my $prefix = q{};
my $help;

my $seqname  = q{};
my $start_nt = q{};
my $end_nt   = q{};

GetOptions ( 'input=s{,}' => \@input_files,
             'prefix=s',  => \$prefix,
             'ascending'  => \$ascending,
             'help'       => \$help,         );

if ( $help or (! @input_files ) ) { 
    die "Format: gff2blocks.pl",
        " --input|-i [input GFF-like files]",
        " --prefix|-p [strip this prefix from all sequence names, e.g., 'chr']",
        " --ascending|-a [optionally, enforce ascending nt order]",
        " --help|-h [this message]",
        "\n",
        ;
}

foreach my $infile (@input_files) { 
    open my $INFILE, '<', $infile or die "Can't open input GFF-like file $infile: $!";
    while (my $input = <$INFILE>) { 
        chomp $input;
        if ( $input =~ /\A (\S+) \t [^\t]* \t [^\t]* \t (\d+) \t (\d+) \t /xms ) { 
            $seqname  = $1;
            $start_nt = $2;
            $end_nt   = $3;
            if ($prefix) { 
                $seqname =~ s/\A$prefix//;
            }
            if ( $ascending and ( $start_nt > $end_nt ) ) { 
                ($start_nt, $end_nt) = ($end_nt, $start_nt);
            }
            print "$seqname\t$start_nt\t$end_nt\n";
        }
    }
    close $INFILE or die "Can't close filehandle to GFF-like file $infile: $!";
}

