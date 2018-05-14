#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $i = 0;

# Store very long file names here:
my $smrtwrap = '/mnt/lustre_scratch_2012/schwarz/work/2015/caenogens/smrtanalysis/install/smrtanalysis_2.3.0.140936/smrtcmds/bin/smrtwrap';
my $pbalign  = '/mnt/lustre_scratch_2012/schwarz/work/2015/caenogens/smrtanalysis/install/smrtanalysis_2.3.0.140936/analysis/bin/pbalign';

my %seen = ();

while (my $input = <>) {
    chomp $input;

    # Ensure that the file both exists and is readable:
    if (! -r $input ) {
        die "Can't read this putative file: $input\n";
    }

    # Enforce very stereotypical input:
    # E.g.:
    # /mnt/lustre_scratch_2012/schwarz/work/2015/caenogens/nigoni/quiver/bax_h5/m150523_165930_42137_c100825762550000001823174211251554_s1_p0.1.bax.h5

    if ( $input =~ /\A (\/mnt\/lustre_scratch_2012\/schwarz\/work
                       \/2015\/caenogens\/nigoni\/quiver\/bax_h5\/.+_s1_p0\.)
                     (\d)
                     \.bax\.h5 
                     \z /xms ) {
        my $stem = $1;
        my $j    = $2;

        # Require a strict order of '_p0.1.bax.h5', '_p0.2.bax.h5', and '_p0.3.bax.h5' for each stem.
        if ( $j == 1 ) {
            if ( exists $seen{$stem} ) { 
                die "Redundant stem detected in: $input\n";
            }
            $seen{$stem} = 1;
            $seen{$j}    = 1;
        }
        elsif ( $j == 2 ) { 
            if (! exists $seen{$stem} ) {
                die "Failed to properly record, earlier, the stem detected in: $input\n";
            }
            $seen{$j}    = 1;
        }
        elsif ( $j == 3 ) {
            if ( ( exists $seen{$stem} ) and ( exists $seen{1} ) and ( exists $seen{2} ) ) {
                %seen = ();  # Reinitialize as empty, to deal with next stanza.

                $i++;
                my $file_no      = sprintf "%02i", $i;
                my $next_file_no = sprintf "%02i", ($i + 1);

                my @k = (1, 2, 3);

                my @h5_files = ();
                foreach my $k (@k) {
                    my $h5_file = $stem . $k . '.bax.h5';
                    if (! -r $h5_file ) {
                        die "Failed to predict readable file ($h5_file) for input: $input\n";
                    }
                    push @h5_files, "$h5_file";
                }

                my $fofn_file = "nigoni_bax.h5_filelist_2015.10.29.$file_no.fofn";
                $fofn_file    = safename($fofn_file);

                my $cmp_h5 = "nigoni_bax.h5_filelist_2015.10.29.$file_no.cmp.h5";
                $cmp_h5    = safename($cmp_h5);

                my $progress = "nigoni_bax.h5_filelist_2015.10.29.$file_no.pbalign_progress.txt";
                $progress    = safename($progress);

                my $start_time = "nigoni_bax.h5_filelist_2015.10.29.$file_no.start_time.txt";
                $start_time    = safename($start_time);

                open my $FOFN, '>', $fofn_file;
                foreach my $h5_file (@h5_files) {
                    print $FOFN "$h5_file\n";
                }
                close $FOFN;

                my $qsub_file = "job_nigoni_mhap_haploid_pbalign_2015.10.29.$file_no.sh";
                $qsub_file    = safename($qsub_file);

                my $next_qsub_file = "job_nigoni_mhap_haploid_pbalign_2015.10.29.$next_file_no.sh";
                $next_qsub_file    = safename($next_qsub_file);

                open my $QSUB, '>', $qsub_file;

                print $QSUB '#!/bin/bash -login', "\n";
                print $QSUB '#PBS -l walltime=012:00:00', "\n";
                print $QSUB '#PBS -l nodes=1:ppn=8', "\n";
                print $QSUB '#PBS -l mem=64gb', "\n";
                print $QSUB "#PBS -N $qsub_file\n";
                print $QSUB '#PBS -q main', "\n";
                print $QSUB '#PBS -M ems394@cornell.edu', "\n";
                print $QSUB '#PBS -m abe', "\n";
                print $QSUB '#PBS -A ged', "\n";
                print $QSUB '#PBS -r n', "\n";
                print $QSUB '#PBS -V', "\n";
                print $QSUB 'cd /mnt/lustre_scratch_2012/schwarz/work/2015/caenogens/nigoni/quiver ;', "\n";
                print $QSUB "echo `date` > $start_time ;\n";
                print $QSUB "$smrtwrap $pbalign --verbose --nproc 8",
                            " --forQuiver $fofn_file",
                            " nigoni_mhap_haploid_2015.10.23.01.scf.fasta",
                            " $cmp_h5",
                            ' 2>',
                            "$progress",
                            " ;\n",
                            ;
                print $QSUB "qsub $next_qsub_file ;\n";
                close $QSUB;
            }
            else {
                die "Failed to observe full three-file set for $input\n";
            }
        }
        else {
            die "Can't identify input as being part of a 1/2/3 file set: $input\n";
        }
    }
    else {
        die "Can't recognize correctly formatted input line: $input\n";
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

