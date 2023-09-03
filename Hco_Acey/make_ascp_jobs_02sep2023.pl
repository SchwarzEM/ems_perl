#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $dirs = q{};
my $dest = q{};

use Cwd;

$dirs = $ARGV[0] if $ARGV[0];
$dest = $ARGV[1] if $ARGV[1];

my $data_ref;

my $start_dir = getcwd;
my $i         = 0;

if ( (! $dirs) or (! $dest) ) {
    die "Format: make_ascp_jobs_02sep2023.pl [list of dirs] [ascp destination] => [one sbatch ascp upload per dir]\n";
}

open my $DEST, '<', $dest;
while ( my $destination = <$DEST> ) {
    chomp $destination;
    # a typical 'destination' for ascp will be a mess of odd characters that don't work well in '/\S/' tests; so, no test here
    if ( exists $data_ref->{'destination'} ) {
        die "Only one destination can be specified\n";
    }
    $data_ref->{'destination'} = $destination;
}
close $DEST;

open my $DIRS, '<', $dirs;
while ( my $input_dir = <$DIRS> ) {
    chomp $input_dir;
    if (! -d $input_dir ) {
        die "Cannot see valid input directory: $input_dir\n";
    }
    my $destination = $data_ref->{'destination'};

    $i++;
    my $j = sprintf "%02u", $i;

    my $job_name = "job_Nippo_sra_ascp_2023.09.03.$j.sh" ;
    $job_name = safename($job_name);

    open my $OUT, '>', $job_name;

    print $OUT '#!/bin/bash', "\n";
    print $OUT '#SBATCH --nodes=1', "\n";
    print $OUT '#SBATCH --partition=RM-shared', "\n";
    print $OUT '#SBATCH --time=048:00:00', "\n";
    print $OUT '#SBATCH --ntasks-per-node=2', "\n";
    print $OUT "#SBATCH --job-name=$job_name\n";
    print $OUT '#SBATCH --mail-type=ALL', "\n";
    print $OUT "cd $start_dir ;\n";
    print $OUT 'module load aspera-connect/3.11.0.5 ;', "\n";
    print $OUT 'ascp -i /jet/home/schwarze/aspera.openssh -QT -l100m -k1 ';
    print $OUT "-d $input_dir ";
    print $OUT "$destination ;\n";
    close $OUT;
}
close $DIRS;

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

