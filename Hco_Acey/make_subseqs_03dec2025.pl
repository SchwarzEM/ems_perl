#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $targets = q{};
$targets    = $ARGV[0] if $ARGV[0];

if (! $targets ) {
    die "Format: make_subseqs_03dec2025.pl [TSV list: target FASTAs + Phobius annots] > [seqkit subseqs trimming SigPs]\n"; 
}

open my $TARGETS, '<', $targets;
while (my $input = <$TARGETS>) {
    chomp $input;
    if ( ( $input =~ /\A (\S+) \t/xms ) and (! -e $1 ) ) {
        my $infile  = $1;
        die "Cannot find putative infile: $infile\n";
    }

    if ( $input =~ /\A (\S+) \t n \d+ [-] \d+ c \d+ \/ (\d+) o /xms ) { 
        my $infile  = $1;
        my $start   = $2;
        my $outfile = $infile;

        $outfile =~ s/\.fa\z//;
        $outfile = "$outfile.nosig.fa";
        $outfile = safename($outfile);

        my $rename = "$infile.w_sig";
        $rename    = safename($rename);

        print "seqkit subseq -R -r $start:-1 $infile > $outfile ;\n";
        print "mv -i $infile $rename ;\n";
        
    }
}
close $TARGETS;

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


