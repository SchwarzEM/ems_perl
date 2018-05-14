#!/usr/bin/perl

# cds2fullnames.pl -- Erich Schwarz <emsch@its.caltech.edu>, 12/5/2012.
# Purpose: given a nametable for WSxxx genes, and a file with CDSes or seq. names starting its lines, change CDS to "WBGene\d+", then print results to STDOUT.  (Note that this output can be pipelined through wbg2fullnames.pl to get full names if desired.)

use strict;
use warnings;
use Getopt::Long;

my %cds2name = ();

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
    die 'Format: wbg2fullnames.pl --names|-n [file format: (WBGene\d+ \s+ pubname \s+ seqname)] --input|-i [1+ files or \'-\' stream] --help|-h [print this message]', "\n", ;
}

open my $NAMES, "<", "$names"  or die "Can't open nametable file $names: $!";
while (my $input = <$NAMES>) { 
    chomp $input;
    if ($input =~ / \A 
                    ( WBGene\d+ ) 
                    \s+ 
                      \S+
                    \s+ 
                    ( \S+ ) 
                  /xms) { 
        my ($wbgene, $seqname);
        ($wbgene, $seqname) = ($1, $2);
        $seqname =~ s/[a-z]\z//;
        $cds2name{$seqname} = $wbgene;
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

        # Convert all CDS names at the start of the line to plain WBGene names, if possible.
        # However, it is important to *recognize* suffixed CDSes and replace them with WBGene synonyms of non-suffixed CDSes.

        if ( $input =~ /\A (\S+) \b/xms ) {
            my $cds_orig = $1;
            my $cds      = $cds_orig;
            $cds =~ s/[a-z]\z//;
            if ( exists $cds2name{$cds} ) { 
                $input =~ s/\A$cds_orig/$cds2name{$cds}/;
            }
        }

        # Print result:
        print "$input\n";
    }
    close $INFILE or die "Can't close filehandle to input file $infile: $!";
}

