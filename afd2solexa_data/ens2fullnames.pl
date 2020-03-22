#!/usr/bin/env perl

# ens2fullnames.pl -- Erich Schwarz <ems394@cornell.edu>, 3/21/2020.
# Purpose: given a nametable for ENSxxx\d genes, throughout input text(s), change "ENSxxx\d" to "ENSxxx\d|name", then print results to STDOUT.

use strict;
use warnings;
use autodie;
use Getopt::Long;

my %ens2name = ();

my @infiles = ();
my $names   = q{};
my $help;

GetOptions(
    "input=s{,}" => \@infiles,
    "names=s"    => \$names,
    "help"       => \$help,
);


# Take first argument as source of file names; do streaming edit of rest.

if ($help or (! $names) or (! @infiles) ) { 
    die 'Format: ens2fullnames.pl --names|-n [file format: (ENSxxx\d \s+ pubname )] --input|-i [1+ files or \'-\' stream] --help|-h [print this message]', "\n", ;
}

open my $NAMES, "<", "$names";
while (my $input = <$NAMES>) { 
    chomp $input;
    if ($input =~ / \A 
                    ( ENS\S*G\d\S+\d ) 
                    \t
                    ( \S+ ) 
                  /xms) { 
        my ($gene_id, $gene_name);
        ($gene_id, $gene_name) = ($1, $2);
        $ens2name{$gene_id} = $gene_id . q{|} . $gene_name ;
    }
    else {
       die "Cannot parse input from $names: $input\n";
    }
}
close $NAMES or die "Can't close filehandle to nametable file $names: $!";

my $INFILE;
foreach my $infile (@infiles) { 
    if ($infile eq '-') {
        # Special case: get the stdin handle
        $INFILE = *STDIN{IO};
    }
    else {
        # Standard case: open the file
        open $INFILE,  '<', $infile or die "Cannot open input file $infile: $!";
    }
    while (my $input = <$INFILE>) { 
        chomp $input;

        # Expand all instances found of WBGene\d+ names in the line; make no effort to check for errors!
        $input =~ s/(ENS\S*G\d\S+\d)\b/$ens2name{$1}/g;

        # Print result:
        print "$input\n";
    }
    close $INFILE or die "Can't close filehandle to input file $infile: $!";
}

