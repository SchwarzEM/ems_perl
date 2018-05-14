#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;

    # sample input line:
    # /mnt/ls15/scratch/users/emsch/work_rsync/2015/adrienne/rsem_2016/replicate_1/atml1_3_rep1.fq	144	60
    #
    # desired output line:
    # /mnt/ls15/scratch/users/emsch/work_rsync/2015/adrienne/rsem_2016/replicate_exact_1/atml1_3_rep1.trim_exact4nt.fq	144	60

    if ( $input =~ /\A (\S+) (\t \d+ \t \d+) \z/xms ) {
        my $file   = $1;
        my $params = $2; 
        if ( $file =~ /rsem_2016\/replicate_(\d+)\/\S+\.fq\z/xms ) {
            my $rep = $1;
            $file =~ s/rsem_2016\/replicate_$rep\//rsem_2016\/replicate_exact_$rep\//;
            $file =~ s/\.fq\z/\.trim_exact4nt\.fq/;
            if (! -e $file) {
                die "Proposed file does not exist: $file\n";
            }
            else {
                print "$file$params\n";
            }
        }
        else {
            die "Cannot rename target file: $file\n"
        }
    }
    else {
        die "Cannot parse input line: $input\n";
    }
}

