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
    die "Format: make_acey_salmon_2022.02.12.02.pl\n",
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

        if ( $input_data_1 !~ / _1 \.(?:fastq|fq) (?:\.gz){0,1} \z/xms ) {
            die "Cannot accept data as a valid FastQ file: $input_data_1\n";
        }

        my $input_data_2 = $input_data_1;
        $input_data_2 =~ s/_1(\.(?:fastq|fq)(?:\.gz){0,1})\z/_2$1/;

        $i++;
        my $j = sprintf "%02u", $i;

        my $output_script = "job_acey_salmon_2022.02.12.$j.sh";
        $output_script    = safename($output_script);

        my $output_prefix = $prefix . ".salmon.2022.02.12.$j";

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
    print $OUT '#SBATCH --time=04:00:00', "\n";
    print $OUT '#SBATCH --ntasks-per-node=16', "\n";
    print $OUT '#SBATCH --job-name=', $output_script, "\n";
    print $OUT '#SBATCH --mail-type=ALL', "\n";

    print $OUT "cd $curr_dir ;\n";

    print $OUT '. $PROJECT/anaconda3/etc/profile.d/conda.sh ;', "\n";
    print $OUT 'conda activate salmon_1.6.0 ;', "\n";

    print $OUT "salmon --no-version-check quant --threads 16 --seqBias --gcBias --posBias --discardOrphansQuasi ";
    print $OUT '--index $PROJECT/Acey/2022.02.13/salmon/dbs/Acey.v2_WBPS16_gentrome_index ';
    print $OUT "--libType A ";
    print $OUT "--mates1 $input_data_1 ";
    print $OUT "--mates2 $input_data_2 ";
    print $OUT "--output $output_prefix ";
    print $OUT '--geneMap $PROJECT/Acey/2022.02.13/annots/Acey_v2.2022.02.15.01.cds2gene.tsv.txt ;';

    print $OUT "\n";

    print $OUT "conda deactivate ;\n";

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

