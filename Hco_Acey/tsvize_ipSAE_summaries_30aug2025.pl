#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Basename;

my $infiles = q{};

$infiles = $ARGV[0] if $ARGV[0];

if (! $infiles ) {
    die "Format: tsvize_ipSAE_summaries_30aug2025.pl [list of 1+ ipSAE summary *.txt files] => 1+ TSV outfiles in working directory\n";
}

open my $INFILES, '<', $infiles;
while ( my $infile = <$INFILES> ) {
    chomp $infile;

    if (! -r $infile ) {
        die "Cannot read putative infile $infile\n";
    }

    my $outfile = basename($infile);
    $outfile =~ s/\.txt\z/.tsv.txt/;

    print "cat $infile";
    print ' | perl -ne \' if( /\S/xms ) { print; } \' | perl -ne \' s/[ ]+/\t/g; print; \' > ';
    print "$outfile ;\n";
}

sub safename {
    my $filename = $_[0];
    my $orig_filename = $filename;
    if (-e $orig_filename) {
        my $suffix1 = 1;
        $filename = $filename . ".$suffix1";
        while (-e $filename) {
            $suffix1++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix1";
        }
    }
    return $filename;
}

