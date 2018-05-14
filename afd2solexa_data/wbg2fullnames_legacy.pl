#!/usr/bin/perl

# wbg2fullnames_legacy.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/20/2008.
# Purpose: change "WBGene\d+" to "WBGene\d+|cgc|seq", and '|' to ';".

use strict;
use warnings;

my %wbg2name = ();

# Take first argument as source of file names; do streaming edit of rest.

if ($#ARGV < 1) { 
    die 'Format: wbg2fullnames_legacy.pl [nametable (WBGene\d+ \s+ pubname \s+ seqname)]  [1+ streams]', "\n", ;
}

my $names = shift @ARGV;

open my $NAMES, "<", "$names" 
    or die "Can't open nametable file $names: $!";

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
            $wbg2name{$wbgene} = '|'. $pubname;
        }
        if ( $pubname ne $seqname ) { 
            $wbg2name{$wbgene} = '|'. $pubname . '|'. $seqname;
        }
    }
}
close $NAMES;

while (my $input = <>) { 
    chomp $input;
    $input =~ s/\|(WBGene\d+)/;$1/g;
    $input =~ s/(WBGene\d+)/$1$wbg2name{$1}/g;
    print "$input\n";
}

