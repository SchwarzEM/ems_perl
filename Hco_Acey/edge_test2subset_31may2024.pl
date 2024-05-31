#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Basename;

my $infile_list = q{};

$infile_list = $ARGV[0] if $ARGV[0];

if (! $infile_list ) {
    die "Format: edge_test2subset_31may2024.pl [list of input edgeR result files] => [subset edgeR files with just Gene, logFC, and FDR]\n"
}

open my $INFILE_LIST, '<', $infile_list;
while (my $infile = <$INFILE_LIST>) {
    chomp $infile;
    if (! -r $infile ) {
        die "Cannot read infile: $infile\n";
    }

    my $basename = basename($infile);
    my $stem     = q{};
    if ( $basename =~ /\A (\S+)_edgeR_exactTest_\d+\.\d+\.\d+\.\d+\.tsv.txt \z/xms ) {
        $stem = $1;
    }
    else {
        die "Cannot parse basename: $basename\n";
    }

    my $outfile = "$stem.sub.tsv.txt";
    $outfile    = safename($outfile);
    open my $OUTFILE, '>', $outfile;

    open my $INFILE, '<', $infile;
    while (my $input = <$INFILE>) {
        chomp $input;
        if ( $input !~ /\A (?: \S+ \t){4} \S+ \z/xms ) {
            die "Mis-tabbed line in input file $infile: $input\n";
        }
        elsif ( $input =~ /\A (\S+) \t (\S+) \t \S+ \t \S+ \t (\S+) \z/xms  ) {
            my $gene  = $1;
            my $logFC = $2;
            my $fdr   = $3;
            print $OUTFILE "$gene\t$logFC\t$fdr\n";
        }
        else {
            die "From input file $infile, cannot parse: $input\n";
        }
    }
    close $INFILE;
    close $OUTFILE;
}
close $INFILE_LIST;

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

