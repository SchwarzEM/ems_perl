#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Basename;
use File::Spec::Functions;

while (my $input_F = <>) {
    chomp $input_F;
    # sample input: fasta_orig/10-1-M13F.fa

    my $dirname    = dirname($input_F);
    my $basename_F = basename($input_F);

    my $basename_R = q{};
    my $tag        = q{};
    if ( $basename_F =~ /\A(\S+)F\.fa/xms ) { 
        $tag        = $1;
        $basename_R = $tag . 'R.fa';
    }
    else {
       	die "Cannot parse F orientation	in $input_F\n";
    }

    my $input_R =  catfile($dirname, $basename_R);

    if (! -r $input_F) {
        die "Cannot read F-input file: $input_F\n";
    }
    if (! -r $input_R) {
       	die "Cannot read R-input file: $input_R\n";
    } 
    my $out_aln = "$tag.merged.aln";
    my $out_seq = "$tag.merged.fa";

    print "merger -asequence $input_F -bsequence $input_R -sreverse2 -outfile $out_aln -outseq $out_seq ;\n";
}

