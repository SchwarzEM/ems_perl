#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Spec::Functions;  # catdir
use Scalar::Util qw(looks_like_number);
use Cwd;

my $i        = 0;
my $workdir  = getcwd;
my $job_stem = "job_renaming_2022.08.01";

while (my $input = <>) {
    # Sample input:
    # N001/i19-j17_1stGATTGTCC-2ndTCGGATTC_Dillman-N01_S575_L002_R1_001.filt1.fq
    chomp $input;
    if ( $input =~ /\A (N\d+) \/ (.+ R1 .+ \.fq) \z/xms ) {
        my $indir    = $1;
        my $infile_1 = $2;

        my $infile_2 = $infile_1;
        $infile_2 =~ s/R1/R2/g;

        my $outfile_1 = $infile_1;
        $outfile_1 =~ s/\.filt1\.fq\z/.filt1.retro.fq/;
       	my $outfile_2 =	$infile_2;
       	$outfile_2 =~ s/\.filt1\.fq\z/.filt1.retro.fq/;

        my $mergefile = $outfile_1;
        $mergefile    =~ s/R1/R1.R2/g;

        my $outdir = "$indir.rev";

        $infile_1 = catfile($indir, $infile_1);        
       	$infile_2 = catfile($indir, $infile_2);

       	$outfile_1 = catfile($outdir, $outfile_1);
        $outfile_2 = catfile($outdir, $outfile_2);
        $mergefile = catfile($outdir, $mergefile);

        $outfile_1 = safename($outfile_1);
        $outfile_2 = safename($outfile_2);
        $mergefile = safename($mergefile);

        if (! -r $infile_1 ) {
            die "Cannot read infile 1: $infile_1\n";
        }
        if (! -r $infile_2 ) {
       	    die	"Cannot read infile 2: $infile_2\n";
       	}

        $i++;
        my $job_no = sprintf "%03i", $i;
        my $job_name = "$job_stem.$job_no.sh";
        $job_name = safename($job_name);

        open my $JOB1, '>', $job_name;        

        print $JOB1 '#!/bin/bash', "\n";
        print $JOB1 '#SBATCH --nodes=1', "\n";
        print $JOB1 '#SBATCH --partition=RM-shared', "\n";
        print $JOB1 '#SBATCH --time=004:00:00', "\n";
        print $JOB1 '#SBATCH --ntasks-per-node=1', "\n";
        print $JOB1 '#SBATCH --job-name=', "$job_name\n";
        print $JOB1 '#SBATCH --mail-type=ALL', "\n";
        print $JOB1 "cd $workdir ;\n";
        print $JOB1 '. $PROJECT/anaconda3/etc/profile.d/conda.sh ;', "\n";
        print $JOB1 "cat $infile_1 | perl -ne ", q{' s/\A([@]\S+) 1/$1\#0\/1 1/; print; '}, " > $outfile_1 ;\n";
        print $JOB1 "cat $infile_2 | perl -ne ", q{' s/\A([@]\S+) 1/$1\#0\/2 2/; print; '}, " > $outfile_2 ;\n";
        print $JOB1 'conda activate seqtk_1.3 ;', "\n";
        print $JOB1 "seqtk mergepe $outfile_1 $outfile_2 > $mergefile ;\n";
        print $JOB1 "conda deactivate ;\n";
        print $JOB1 'conda activate seqkit_2.1.0 ;', "\n";
        print $JOB1 "seqkit stats --all --threads 1 $mergefile > $mergefile.seqkit_stats.txt ;\n";
        print $JOB1 "conda deactivate ;\n";
        close $JOB1;
    }
    else {
        die "Cannot parse input: $input\n";
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

