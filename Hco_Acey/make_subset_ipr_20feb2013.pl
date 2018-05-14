#!/usr/bin/env perl

use strict;
use warnings;
use Cwd;
use File::Basename;
use File::Spec;

my ($day, $month, $year) = (localtime)[3,4,5];

$year  = sprintf "%04d", ($year+1900);
$month = sprintf "%02d", ($month+1);
$day   = sprintf "%02d", $day;

my $current_date = "$year.$month.$day";
my $current_dir  = cwd ;

while (my $input = <>) { 
    chomp $input;
    my $stem = $input;
    $stem    = basename $stem;
    $stem    =~ s/\.fa\z//;

    my $target_dir = $stem . '.dir';
    $target_dir    = File::Spec->catdir( $current_dir, $target_dir );

    my $qsub_script = 'job_iprscan_' . $stem . '_' . $current_date . '.sh';
    system "mkdir $target_dir";
    system "mv $input $target_dir";
    open my $SCRIPT, '>', $qsub_script or die "Can't open qsub script: $qsub_script\n";

    print $SCRIPT '#!/bin/bash -login', "\n";
    print $SCRIPT '#PBS -l walltime=024:00:00', "\n";
    print $SCRIPT '#PBS -l nodes=1:ppn=4', "\n";
    print $SCRIPT '#PBS -l mem=10gb', "\n";
    print $SCRIPT "#PBS -N $qsub_script\n";
    print $SCRIPT '#PBS -q main', "\n";
    print $SCRIPT '#PBS -M ems394@cornell.edu', "\n";
    print $SCRIPT '#PBS -m abe', "\n";
    print $SCRIPT '#PBS -A ged-intel11', "\n";
    print $SCRIPT '#PBS -r n', "\n";
    print $SCRIPT '#PBS -V', "\n";
    print $SCRIPT "cd $target_dir ;\n";
    print $SCRIPT 'module load iprscan/4.8 ;', "\n";
    print $SCRIPT 'iprscan -cli -email ems394@cornell.edu', 
                  " -i $input -o $stem.iprscan.raw -format raw -iprlookup -goterms ;\n", 
                  ;

    close $SCRIPT or die "Can't close filehandle to qsub script: $qsub_script\n";
    system "mv $qsub_script $target_dir";
} 

