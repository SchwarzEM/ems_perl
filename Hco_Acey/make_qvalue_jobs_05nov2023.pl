#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Basename;

while (my $input = <>) {
    chomp $input;
    my $output = basename($input);
    my $error  = basename($input);
    $output =~ s/\.tsv\.txt\z/.stats.raw.txt/;
    $error  =~ s/\.tsv\.txt\z/.stats.err/;

    $output = safename($output);
    $error  = safename($error);

    print '$PROJECT/ems_perl/Hco_Acey/motif_group_fisher_03nov2023.pl ',
          "$input ",
          '$PROJECT/mambaforge-pypy3/envs/meme_5.4.1/bin/qvalue ',
          "1>$output 2>$error ;\n",
          ;
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

