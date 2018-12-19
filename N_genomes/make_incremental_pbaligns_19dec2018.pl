#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Basename;

my $i = 0;

# Store very long file names here:
my $smrtwrap = '/pylon5/mc5fpnp/schwarze/smrtanalysis_precursors/install/smrtanalysis_2.3.0.140936/smrtcmds/bin/smrtwrap';
my $pbalign  = '/pylon5/mc5fpnp/schwarze/smrtanalysis_precursors/install/smrtanalysis_2.3.0.140936/analysis/bin/pbalign';

my $suffix          = q{};
my $genome_sequence = q{};
my $working_dir     = q{};

my %seen = ();

while (my $input = <@ARGV>) {
    chomp $input;

    # Ensure that any file being invoked both exists and is readable:
    if ( ($suffix) and (! -r $input ) ) {
        die "Can't read this putative file: $input\n";
    }

    # Accept the very first file as a genome sequence, 
    # but require all further files to be '_p0.1.bax.h5', '_p0.2.bax.h5', and '_p0.3.bax.h5' for each stem.
    # Enforce fairly stereotypical input:
    # E.g.:
    # /mnt/lustre_scratch_2012/schwarz/work/2015/caenogens/nigoni/quiver/bax_h5/m150523_165930_42137_c100825762550000001823174211251554_s1_p0.1.bax.h5

    # Allow a specific suffix to be supplied, ONCE, making it easier to tag distinct parallel jobs.
    if (! $suffix ) {
        $suffix = $input;
        if ( $suffix !~ /\A \S+ \z/xms ) {
            die "Suffix (\"$suffix\") cannot have space characters\n";
        }
    }

    # Allow a genome sequence to be supplied, ONCE; require it be a full filename.
    elsif (! $genome_sequence ) {
        $genome_sequence = $input;

        # Use this to extract the working directory for the qsub-able scripts, while keeping only the filename for the genome sequence file.
        $working_dir     = dirname($genome_sequence);
        $genome_sequence = basename($genome_sequence);
    }

# /mnt/ls15/scratch/users/emsch/work_rsync/2015/caenogens/tropicalis/quiver/tropicalis/*/*/Analysis_Results/*.bax.h5
# /pylon5/mc5fpnp/schwarze/wallacei/tropicalis_data/tropicalis/*/*/Analysis_Results/*.bax.h5

    elsif ( $input =~ /\A ( \/pylon5\/mc5fpnp\/schwarze\/wallacei\/tropicalis_data \/ \S+ \/ [^\/\s]+ _s1_p0\.)
                     (\d)
                     \.bax\.h5 
                     \z /xms ) {
        my $stem = $1;
        my $j    = $2;

        # enforce pre-existence of a define genome sequence file:
        if (! $genome_sequence ) {
            die "Need to list full genome sequence (with full pathdir) before listing *_p0.[1-3].bax.h5 files\n";
        }

        # Require a strict order of '_p0.1.bax.h5', '_p0.2.bax.h5', and '_p0.3.bax.h5' for each stem.
        if ( $j == 1 ) {
            if ( exists $seen{$stem} ) { 
                die "Redundant stem detected in: $input\n";
            }
            $seen{$stem} = 1;
            $seen{$j}    = 1;
        }
        elsif ( $j == 2 ) { 
            if ( (! exists $seen{$stem} ) or (! exists $seen{1} ) ) {
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

                my @k_vals = (1, 2, 3);

                my @h5_files = ();
                foreach my $k (@k_vals) {
                    my $h5_file = $stem . $k . '.bax.h5';
                    if (! -r $h5_file ) {
                        die "Failed to predict readable file ($h5_file) for input: $input\n";
                    }
                    push @h5_files, "$h5_file";
                }

                my $fofn_file = "tropicalis.bax.h5_filelist_" . "$suffix.$file_no.fofn";
                $fofn_file    = safename($fofn_file);

                my $cmp_h5 = "tropicalis.bax.h5_filelist_" . "$suffix.$file_no.cmp.h5";
                $cmp_h5    = safename($cmp_h5);

                my $progress = "tropicalis.bax.h5_filelist_" . "$suffix.$file_no.pbalign_progress.txt";
                $progress    = safename($progress);

                my $start_time = "tropicalis.bax.h5_filelist_" . "$suffix.$file_no.start_time.txt";
                $start_time    = safename($start_time);

                open my $FOFN, '>', $fofn_file;
                foreach my $h5_file (@h5_files) {
                    print $FOFN "$h5_file\n";
                }
                close $FOFN;

                my $qsub_file = "job_tropicalis.pbalign_" . "$suffix.$file_no.sh";
                $qsub_file    = safename($qsub_file);

                my $next_qsub_file = "job_tropicalis.pbalign_" . "$suffix.$next_file_no.sh";
                $next_qsub_file    = safename($next_qsub_file);

                open my $QSUB, '>', $qsub_file;

                print $QSUB '#!/bin/bash', "\n";
                print $QSUB '#SBATCH --nodes=1', "\n";
                print $QSUB '#SBATCH --partition=RM-shared', "\n";
                print $QSUB '#SBATCH --time=012:00:00', "\n";
                print $QSUB '#SBATCH --ntasks-per-node=8', "\n";
                print $QSUB '#SBATCH --constraint=EGRESS', "\n";
                print $QSUB "#SBATCH --job-name=$qsub_file\n";
                print $QSUB '#SBATCH --mail-type=ALL', "\n";
                print $QSUB "cd $working_dir ;\n";
                print $QSUB 'module load gcc/4.8.4 ;', "\n";
                print $QSUB "echo `date` > $start_time ;\n";
                print $QSUB "$smrtwrap $pbalign --verbose --nproc 8",
                            " --forQuiver $fofn_file",
                            " $genome_sequence",
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

