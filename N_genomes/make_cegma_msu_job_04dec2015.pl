#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Cwd;
use File::Spec::Functions;

my $basedir = getcwd;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \.fa \z/xms ) {
        my $stem = $1;

        my $dir = $stem . '_dir';
        $dir    = catfile($basedir, $dir);
        my $tmp = catfile($dir, 'tmp');

        # To make Perl's weird syntax work here, there needs to be a leading negation ('!')
        !system "mkdir -p $tmp" or die "Cannot mkdir \"mkdir -p $tmp\"\n";
        !system "mv $input $dir" or die "Cannot \"mv $input $dir\"\n";

        my $cegma = $stem . '_cegma';

        my $script = 'job_cegma_' . $stem . '_2015.12.04.sh';
        $script    = safename($script);

        open my $SCRIPT, '>', $script;

        print $SCRIPT '#!/bin/bash -login', "\n";
        print $SCRIPT '#PBS -l walltime=006:00:00', "\n";
        print $SCRIPT '#PBS -l nodes=1:ppn=8', "\n";
        print $SCRIPT '#PBS -l mem=32gb', "\n";
        print $SCRIPT "#PBS -N $script\n";
        print $SCRIPT '#PBS -q main', "\n";
        print $SCRIPT '#PBS -M ems394@cornell.edu', "\n";
        print $SCRIPT '#PBS -m abe', "\n";
        print $SCRIPT '#PBS -A ged', "\n";
        print $SCRIPT '#PBS -r n', "\n";
        print $SCRIPT '#PBS -V', "\n";
        print $SCRIPT "cd $dir ;\n";
        print $SCRIPT 'module load CEGMA/2.4 ;', "\n";
        print $SCRIPT 'export CEGMATMP=tmp ;', "\n";
        print $SCRIPT "cegma --genome $input --output $cegma --threads 8 ;\n";

        close $SCRIPT;
    }
    else {
        die "Can't parse input: $input\n";
    }
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

