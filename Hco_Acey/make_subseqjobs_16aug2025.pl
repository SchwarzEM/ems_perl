#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $file2targ  = q{};
my $targ2range = q{};

$file2targ  = $ARGV[0] if $ARGV[0];
$targ2range = $ARGV[1] if $ARGV[1];

my $data_ref;
my $header = 1;

if ( (! $file2targ ) or (! $targ2range ) ) {
    die "Format: make_subseqjobs_16aug2025.pl [file to target TSV] [target to range TSV] > [batch job of seqkit subseqs]\n"

}

open my $FILE2TARG, '<', $file2targ;
while (my $input = <$FILE2TARG>) {
    chomp $input;
    # Sample input: 
    # seqs/CD122.pep.fa       CD122
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
        my $file   = $1;
        my $target = $2;
        if ( exists $data_ref->{'target'}->{$target} ) {
            die "For file2targ $file2targ, redundant target annotations at: $input\n";
        }
        $data_ref->{'target'}->{$target}->{'file'} = $file;
    }
    else {
        die "From file2targ $file2targ, cannot parse: $input\n";
    }
}
close $FILE2TARG;

open my $TARG2RANGE, '<', $targ2range;
while (my $input = <$TARG2RANGE>) {
    chomp $input;
    # Sample input: 
    # CD122   27-240
    if ( $input =~ /\A (\S+) \t (\d+) [-] (\d+) \z/xms ) {
        my $target = $1;
        my $start  = $2;
        my $stop   = $3;
        if ( ( exists $data_ref->{'target'}->{$target}->{'start'} ) or ( exists $data_ref->{'target'}->{$target}->{'stop'} ) ) {
            die "For targ2range $targ2range, redundant target annotations at: $input\n";
        }
        $data_ref->{'target'}->{$target}->{'start'} = $start;
        $data_ref->{'target'}->{$target}->{'stop'}  = $stop;
    }
    else {
        die "From file2targ $file2targ, cannot parse: $input\n";
    }
}
close $TARG2RANGE;

my @targets = sort keys %{ $data_ref->{'target'}  };
foreach my $target (@targets) {
    my $infile = $data_ref->{'target'}->{$target}->{'file'};
    my $start  = $data_ref->{'target'}->{$target}->{'start'};
    my $stop   = $data_ref->{'target'}->{$target}->{'stop'};

    my $outfile = $infile;
    my $suffix  = '_' . $start . '-' . $stop;
    $outfile =~ s/\.pep/$suffix.pep/;

    if ($header) {
        print '#!/bin/bash', "\n";
        print '#SBATCH --nodes=1', "\n";
        print '#SBATCH --partition=RM-shared', "\n";
        print '#SBATCH --time=024:00:00', "\n";
        print '#SBATCH --ntasks-per-node=2', "\n";
        print '#SBATCH --job-name=job_Heligo_subseq_2025.08.16.01.sh', "\n";
        print '#SBATCH --mail-type=ALL', "\n";
        print 'cd $PROJECT/heligmo ;', "\n";
        print 'source $HOME/.bashrc_mamba ;', "\n";
        print '. $PROJECT/mambaforge-pypy3/etc/profile.d/mamba.sh ;', "\n";
        print 'mamba activate seqkit_2.7.0 ;', "\n";
        $header = 0;
    }
    print "seqkit subseq -R -r $start:$stop $infile > $outfile ;\n"
}

print "mamba deactivate ;\n";
