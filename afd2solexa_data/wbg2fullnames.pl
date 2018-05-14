#!/usr/bin/perl

# wbg2fullnames.pl -- Erich Schwarz <emsch@its.caltech.edu>, 12/5/2012.
# Purpose: given a nametable for WSxxx genes, throughout input text(s), change "WBGene\d+" to "WBGene\d+|seq|cgc", then print results to STDOUT.

use strict;
use warnings;
use Getopt::Long;

my %wbg2name = ();

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
                    ( \S+ ) 
                    \s+ 
                    ( \S+ ) 
                  /xms) { 
        my ($wbgene, $pubname, $seqname);
        ($wbgene, $pubname, $seqname) = ($1, $2, $3);
        $seqname =~ s/[a-z]\z//;
        $wbg2name{$wbgene} = $seqname;
        if ( $pubname eq $seqname ) {
            $wbg2name{$wbgene} = $wbgene . q{|} . $pubname;
        }
        if ( $pubname ne $seqname ) { 
            $wbg2name{$wbgene} = $wbgene . q{|} . $seqname . q{|} . $pubname;
        }
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
        $input =~ s/(WBGene\d+)/$wbg2name{$1}/g;

        # Print result:
        print "$input\n";
    }
    close $INFILE or die "Can't close filehandle to input file $infile: $!";
}

