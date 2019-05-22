#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Cwd;
use File::Spec::Functions;  # catdir, catfile

my $start_dir = getcwd;
my $lib_no    = q{};
my $init_desc = q{};
my $lib_desc  = q{};
my $recheck1  = 0;

my $header      = 1;
my $output_file = q{};

while (my $input = <>) {
    chomp $input;

    # Sample three input lines:
    # 
    # Lane #1 : (17751) Index #13 Acey intestines 19dpi -DEX biorep 1
    #      https://jumpgate.caltech.edu/library/17751
    # https://jumpgate.caltech.edu/runfolders/volvox02/161005_SN787_0580_BH3WGFBCXY/Unaligned/Project_17751_index13/Sample_17751/

    if ( (! $lib_no) and (! $init_desc) and (! $lib_desc) and (! $recheck1) 
          and ( $input =~ /\A Lane [ ] [#]\d+ [ ] : [ ] \( (\d+) \) [ ] Index [ ] [#]\d+ [ ] (\S.*\S) \s* \z/xms ) ) { 
        $lib_no    = $1;
        $init_desc = $2;

        $lib_desc = $init_desc;
        $lib_desc =~ s/\s+/_/g;
        $lib_desc =~ s/non[-]intestines/non.intestines/g;
        $lib_desc =~ s/[-]/no/g;
    }
    elsif ( $lib_no and $init_desc and $lib_desc and (! $recheck1) 
            and ( $input =~ /\A [ ]+ https[:]\/\/jumpgate\.caltech\.edu\/library\/ (\d+) \s* \z/xms ) ) { 
        my $putative_lib_no1 = $1;
        if ( $putative_lib_no1 != $lib_no ) {
            die "Inconsistent library number in: $input\n";
        }
        $recheck1 = 1;
    }
    elsif ( $lib_no and $init_desc and $lib_desc and $recheck1 
            and ( $input =~ /\A https:\/\/jumpgate\.caltech\.edu
                            \/runfolders\/volvox02
                            \/ [a-zA-Z0-9]+ _ [a-zA-Z0-9]+ _ [a-zA-Z0-9]+ _ [a-zA-Z0-9]+
                            \/Unaligned\/Project_ (\d+) _index \d+ \/Sample_ (\d+) \/ \z/xms ) ) {
        my $putative_lib_no1 = $1;
        my $putative_lib_no2 = $2;
        if ( ( $putative_lib_no1 != $lib_no ) or ( $putative_lib_no2 != $lib_no ) ) {
            die "Inconsistent library number in: $input\n";
        }

        # If we get through all of these, we can generate a stanza of commands.
        my $target_subdir = $lib_no . q{_} . $lib_desc;
        $target_subdir = catdir($start_dir, $target_subdir);
        if (! -r $target_subdir ) {
            die "Cannot read target subdirectory: $target_subdir\n";
        }

        my $target_raw_reads = "$target_subdir" . q{/} . "$lib_no" . '_*.fastq.gz';

        my $target_output1 = "$target_subdir.filt1.fq";
        $target_output1    = safename($target_output1);

        # print header lines only once, at the top of the output script

        print '#!/bin/bash', "\n" if $header;
        print '#SBATCH --nodes=1', "\n" if $header;
        print '#SBATCH --partition=RM-shared', "\n" if $header;
        print '#SBATCH --time=48:00:00', "\n" if $header;
        print '#SBATCH --ntasks-per-node=1', "\n" if $header;
        print '#SBATCH --constraint=EGRESS', "\n" if $header;
        print '#SBATCH --job-name=job_Acey_rnaseq_2019.05.22.01.sh', "\n" if $header;
        print '#SBATCH --mail-type=ALL', "\n" if $header;
        print "cd $start_dir ;\n" if $header;

        $header = 0;

        print "zcat $target_raw_reads | quality_trim_fastq.pl -q 33 -u 50 -i - -o $target_output1 ;\n";

        # zero out these values:
        $lib_no    = q{};
        $init_desc = q{};
        $lib_desc  = q{};
        $recheck1  = 0;
    }
    else { 
        die "Cannot parse input line: $input\n";
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

