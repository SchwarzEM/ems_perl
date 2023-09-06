#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $infile = q{};
$infile    = $ARGV[0] if $ARGV[0];

if (! $infile ) {
    die "Format: make_biosample_tablings_05sep2023.pl [input directory list] > [listing line commands]\n";
}

open my $INFILE, '<', $infile;
while ( my $input = <$INFILE> ) {
    chomp $input;
    if (! -d $input ) {
        die "Not detectable as a real directory: $input\n";
    }

    my $output = "$input.annots.orig.tsv.txt";
    $output    = safename($output);

    print 'tab_biosample.reads_2023.09.03.02.pl <(ls ';
    print "$input";
    print '/*) > ';
    print "$output ;\n";
}
close $INFILE;

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

