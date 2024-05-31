#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Getopt::Long;
use Cwd;
use List::MoreUtils qw(uniq);

my $infile      = q{};
my $curr_dir    = getcwd;
my $i           = 1;
my %file2script = ();

my @prefixes = ();

my $data_ref;

my $help;

GetOptions ( 'infiles=s' => \$infile,
             'start=i'   => \$i,
             'help'      => \$help,   );

if ( $help or (! $infile) ) {
    die "Format: make_nippo_salmon_2024.05.30.01.pl\n",
        "    --infile|-i   <input table, with prefix-file(s) pairs>\n",
        "    --start|-s    <starting index; default of 1>\n",
        "    --help|-h     [print this message]\n",
        ;
}

$i--;

open my $INFILE, '<', $infile;

while (my $input = <$INFILE>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) { 
        my $prefix       = $1;
        my $input_data_1 = $2;

        # .R1.filt1.fq
        if ( $input_data_1 !~ / R1 \.filt1 \.(?:fastq|fq) (?:\.gz){0,1} \z/xms ) {
            die "Cannot accept data as a valid FastQ file: $input_data_1\n";
        }

        my $input_data_2 = $input_data_1;
        $input_data_2 =~ s/R1(\.filt1\.(?:fastq|fq)(?:\.gz){0,1})\z/R2$1/;

        if (! -e $input_data_1 ) {
            die "R1 input file does not exist: $input_data_1\n";
        }

        if (! -e $input_data_2 ) {
            die "R2 input file does not exist: $input_data_2\n";
        }

        $i++;
        my $j = sprintf "%03u", $i;

        my $output_script = "job_nippo_salmon_2024.05.30.$j.sh";
        $output_script    = safename($output_script);

        my $output_prefix = $prefix . ".salmon_2024.05.30.$j";

        if ( exists $data_ref->{'prefix'}->{$prefix} ) {
            die "Redundant input prefix: $prefix\n";
        }

        push @prefixes, $prefix;

        $data_ref->{'prefix'}->{$prefix}->{'input_data_1'}  = $input_data_1;
        $data_ref->{'prefix'}->{$prefix}->{'input_data_2'}  = $input_data_2;
        $data_ref->{'prefix'}->{$prefix}->{'output_script'} = $output_script; 
        $data_ref->{'prefix'}->{$prefix}->{'output_prefix'} = $output_prefix;
    }
}    

$i = 1;

foreach my $prefix (@prefixes) {
    my $input_data_1  = $data_ref->{'prefix'}->{$prefix}->{'input_data_1'};
    my $input_data_2  = $data_ref->{'prefix'}->{$prefix}->{'input_data_2'};
    my $output_script = $data_ref->{'prefix'}->{$prefix}->{'output_script'};
    my $output_prefix = $data_ref->{'prefix'}->{$prefix}->{'output_prefix'};

    open my $OUT, '>', $output_script;

    print $OUT '#!/bin/bash', "\n";
    print $OUT '#SBATCH --nodes=1', "\n";
    print $OUT '#SBATCH --partition=RM-shared', "\n";
    print $OUT '#SBATCH --time=48:00:00', "\n";
    print $OUT '#SBATCH --ntasks-per-node=16', "\n";
    print $OUT '#SBATCH --job-name=', $output_script, "\n";
    print $OUT '#SBATCH --mail-type=ALL', "\n";

    print $OUT "cd $curr_dir ;\n";
    print $OUT 'source $HOME/.bashrc_mamba ;', "\n";
    print $OUT '. $PROJECT/mambaforge-pypy3/etc/profile.d/mamba.sh ;', "\n";
    print $OUT 'mamba activate salmon_1.10.2 ;', "\n";

    print $OUT "salmon --no-version-check quant --numBootstraps 100 --threads 16 --seqBias --gcBias --posBias --discardOrphansQuasi ";
    print $OUT '--index $PROJECT/nippo_steiner/Nippo_genome/2024.02.19/salmon/dbs/Nippo_2023.07.15.01_gentrome_index ';
    print $OUT "--libType A ";
    print $OUT "--mates1 $input_data_1 ";
    print $OUT "--mates2 $input_data_2 ";
    print $OUT "--output $output_prefix ";
    print $OUT '--geneMap $PROJECT/nippo_steiner/Nippo_genome/2024.02.19/annots/nippo_braker_combined_2023.07.15.cds2gene.tsv.txt ;';

    print $OUT "\n";

    print $OUT "mamba deactivate ;\n";

    if ( exists $prefixes[$i] ) {
        my $next_output_script = $data_ref->{'prefix'}->{$prefixes[$i]}->{'output_script'};
        print $OUT "sbatch $next_output_script ;\n";
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

