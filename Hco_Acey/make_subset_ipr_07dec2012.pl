#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my @infiles    = ();
my $target_dir = q{};
my $sublabel   = q{};
my $help;

GetOptions ( 'infiles=s{,}' => \@infiles,
             'target=s'     => \$target_dir,
             'sublabel=s'   => \$sublabel,
             'help'         => \$help,     );

if ( $help or (! @infiles) or (! $target_dir) or (! $sublabel) ) { 
    die "Format: make_subset_ipr_07dec2012.pl\n",
        "    --infile|-i    <input stream/files>\n",
        "    --target|-t    [target main directory, e.g., /mnt/scratch/emsch/blast2go/iprscan]\n",
        "    --sublabel|-s  [sublabel for jobs/etc., e.g., '16nov2012']\n",
        "    --help|-h      [print this message]\n",
        ;
}

foreach my $input (@infiles) { 
    my $stem = $input;
    $stem    =~ s/\.fa\z//;

    my $sub_target_dir = $stem . '.dir';
    my $qsub_script = 'job_iprscan_' . $stem . q{_} . "$sublabel.sh";
    system "mkdir $sub_target_dir";
    system "mv $input $sub_target_dir";
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
    print $SCRIPT "cd $target_dir/$sub_target_dir ;\n";
    print $SCRIPT 'module load iprscan ;', "\n";
    print $SCRIPT 'iprscan -cli -email ems394@cornell.edu', 
                  " -i $input -o $stem.iprscan.raw -format raw -iprlookup -goterms ;\n", 
                  ;

    close $SCRIPT or die "Can't close filehandle to qsub script: $qsub_script\n";
    system "mv $qsub_script $sub_target_dir";
} 

