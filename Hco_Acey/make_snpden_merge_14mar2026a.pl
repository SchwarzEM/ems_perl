#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my @chrs = qw ( chrI chrII chrIII chrIV chrV chrX );

while ( my $input = <> ) {
    chomp $input;
    if ( $input =~ /\A (\S+) chrI (\S+\.snpden) \z/xms ) {
        my $text1 = $1;
        my $text2 = $2;

        my $outfile = $text1 . 'all_chrs' . $text2;
        $outfile    = safename($outfile);

        my @infiles = ();
        foreach my $chr (@chrs) {
            my $infile = "$text1$chr$text2";
            if (! -e $infile) {
                die "Can't find infile: $infile\n";
            }
            push @infiles, $infile;
        }
        print "cat @infiles ", '| egrep -v "^CHROM" ', "> $outfile ;\n";
    }
    else {
        die "Cannot parse input: $input\n";
    }
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

