#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

my $i = 0;

my @input_archives = ();
my %archive2script = ();

while (my $input = <>) {
    chomp $input;
    if ( $input !~ /\A \S+ \.tar\.gz \z/xms ) { 
        die "Cannot accept input file as a valid archive: $input\n";
    }
    push @input_archives, $input;
}

@input_archives = sort @input_archives;
@input_archives = uniq @input_archives;

foreach my $input_archive (@input_archives) {
    $i++;
    my $j = sprintf "%02u", $i;

    my $output_script = "job_unpack_AFire_new_PacBio_2015.09.30.$j.sh";
    $output_script    = safename($output_script);

    $archive2script{$input_archive} = $output_script;
}

$i = 0;    

foreach my $input_archive (@input_archives) {
    my $k = ($i + 1);

    my $output_script = $archive2script{$input_archive};

    open my $OUT, '>', $output_script;

    print $OUT '#!/bin/bash -login', "\n";
    print $OUT '#PBS -l walltime=024:00:00', "\n";
    print $OUT '#PBS -l nodes=1:ppn=1', "\n";
    print $OUT '#PBS -l mem=32gb', "\n";
    print $OUT "#PBS -N $output_script\n";
    print $OUT '#PBS -q main', "\n";
    print $OUT '#PBS -M ems394@cornell.edu', "\n";
    print $OUT '#PBS -m abe', "\n";
    print $OUT '#PBS -A ged', "\n";
    print $OUT '#PBS -r n', "\n";
    print $OUT '#PBS -V', "\n";
    print $OUT "cd /mnt/ls15/scratch/users/emsch/VC2010/data_25sep2016/pacbio_aug2016 ;\n";
    print $OUT "zcat $input_archive | tar xf - ;\n";

    if ( exists $input_archives[$k] ) {
        my $next_output_archive = $input_archives[$k];
        my $next_output_script  = $archive2script{$next_output_archive};
        print $OUT "qsub $next_output_script ;\n";
    }

    close $OUT;
    $i++;
}

sub safename {
    my $_filename = $_[0];
    my $_orig_filename = $_filename;
    if (-e $_orig_filename) {
        my $_suffix1 = 1;
        $_filename = $_filename . ".$_suffix1";
        while (-e $_filename) {
            $_suffix1++;
            $_filename =~ s/\.\d+\z//xms;
            $_filename = $_filename . ".$_suffix1";
        }
    }
    return $_filename;
}

