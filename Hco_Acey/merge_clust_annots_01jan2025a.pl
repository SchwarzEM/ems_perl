#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Basename;
use List::MoreUtils qw(uniq);

my @infiles = ();

while (my $infile = <>) {
    chomp $infile;
    if (! -r $infile ) {
        die "Cannot read putative input file: $infile\n";
    }
    push @infiles, $infile;
}

@infiles = uniq(@infiles);

foreach my $infile (@infiles) {
    my $basename = basename($infile);
    my $clust_id = q{};
    if ( $basename =~ /\A \S+ \.vs\. (\S+) \.stats\.filt1\.txt \z/xms ) {
        $clust_id = $1;
    }
    else {
        die "Cannot parse basename ($basename) of input file: $infile\n";
    }
    open my $INFILE, '<', $infile;
    while (my $input = <$INFILE>) {
        chomp $input;
        if ( $input !~ /\A Motif /xms ) {
             print "$clust_id\t$input\n";
        }
    }
    close $INFILE;
}

