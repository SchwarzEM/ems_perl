#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my @read_nos      = ();
my @command_lines = ();

my $start_read_no = q{};
my $end_read_no   = q{};

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A raw_reads\/CS004_NoIndex_L002_R1_(\d+)\.fastq.gz \z/xms ) { 
        my $read_no  = $1;
        push @read_nos, $read_no;

        my $r1_reads = $input;
        my $r2_reads = $r1_reads;
        $r2_reads =~ s/CS004_NoIndex_L002_R1_/CS004_NoIndex_L002_R2_/;
        if ( (! -r $r1_reads ) or (! -r $r2_reads ) ) {
            die "Can't read one or both readfiles: $r1_reads or $r2_reads\n";
        }

        my $command1 = 'bowtie2 -p 8 -x nigoni_mhap_pbjelly_quiver.decont_2015.11.11'
                       . " -1 $r1_reads -2 $r2_reads"
                       . " -S meyer_illumina_PE.$read_no.vs.nigoni_mhap_pbjelly_quiver.decont_2015.11.11.sam ;"
                       ;

        my $command2 = "samtools view -b -o meyer_illumina_PE.$read_no.vs.nigoni_mhap_pbjelly_quiver.decont_2015.11.11.bam"
                       . " meyer_illumina_PE.$read_no.vs.nigoni_mhap_pbjelly_quiver.decont_2015.11.11.sam ;"
                       ;

        my $command3 = "mv -i meyer_illumina_PE.$read_no.vs.nigoni_mhap_pbjelly_quiver.decont_2015.11.11.sam sam_files ;";

        my $command4 = "samtools sort -@ 8 -T tmp.$read_no"
                       . " -o meyer_illumina_PE.$read_no.vs.nigoni_mhap_pbjelly_quiver.decont_2015.11.11.sorted.bam"
                       . " meyer_illumina_PE.$read_no.vs.nigoni_mhap_pbjelly_quiver.decont_2015.11.11.bam ;"
                       ;

        my $command5 = "mv -i meyer_illumina_PE.$read_no.vs.nigoni_mhap_pbjelly_quiver.decont_2015.11.11.bam unsorted_bam_files ;";

        my $command6 = "samtools index meyer_illumina_PE.$read_no.vs.nigoni_mhap_pbjelly_quiver.decont_2015.11.11.sorted.bam ;";

        push @command_lines, $command1, $command2, $command3, $command4, $command5, $command6;
    }
    else {
        die "Can't parse input line: $input\n";
    }
}

if (@command_lines) {
    $start_read_no = $read_nos[0];
    $end_read_no   = $read_nos[-1];

    my $start = $start_read_no ;
    $start    =~ s/\A[0]+// ;
    my $end   = $end_read_no;
    $end      =~ s/\A[0]+// ;

    my $interval      = q{};
    if ( $start == $end ) {
        $interval = $start_read_no;
    }
    else { 
        $interval = $start_read_no . q{-} . $end_read_no;
    }

    my $outscript     = 'job_bowtie_nigoni_reads_' . $interval . '_2015.11.14.01.sh';
    $outscript        = safename($outscript);

    open my $OUT, '>', $outscript;

    print $OUT '#!/bin/bash -login', "\n";
    print $OUT '#PBS -l walltime=024:00:00', "\n";
    print $OUT '#PBS -l nodes=1:ppn=8', "\n";
    print $OUT '#PBS -l mem=16gb', "\n";
    print $OUT '#PBS -N', " $outscript\n";
    print $OUT '#PBS -q main', "\n";
    print $OUT '#PBS -M ems394@cornell.edu', "\n";
    print $OUT '#PBS -m abe', "\n";
    print $OUT '#PBS -A ged', "\n";
    print $OUT '#PBS -r n', "\n";
    print $OUT '#PBS -V', "\n";
    print $OUT 'cd /mnt/ls15/scratch/users/emsch/work_rsync/2015/caenogens/nigoni/pilon ;', "\n";
    print $OUT 'module load bowtie2/2.2.6 ;', "\n";
    print $OUT 'module load SAMTools/1.2 ;', "\n";
    foreach my $command (@command_lines) {
        print $OUT "$command\n";
    }
    close $OUT;
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

