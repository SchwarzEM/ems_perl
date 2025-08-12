#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $infile  = q{};
my $outfile = q{};
my $i       = 0;
my $j       = sprintf "%02i", $i;

$infile = $ARGV[0] if $ARGV[0];

if (! $infile ) {
    die "Format split_mailboxes_11aug2025.pl [infile] => [1+ split outfiles]\n";
}

my $INFILE;
my $OUTFILE;

$outfile = "$infile.split.$i.txt";

open $INFILE, '<', $infile;
while ( my $input = <$INFILE> ) {
    chomp $input;
    if ( $input =~ /\A Return-Path:[ ][<] /xms ) {
        $i++;
        $j       = sprintf "%02i", $i;
        $outfile = "$infile.split.$j.txt";
        $outfile = safename($outfile);
        open $OUTFILE, '>', $outfile;
    }
    print $OUTFILE "$input\n";
}
close $INFILE;
close $OUTFILE;

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
