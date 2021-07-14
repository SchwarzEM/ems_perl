#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Cwd;

my $i        = 0;
my $list     = q{};
my $job      = q{};

my $work_dir = getcwd;

$list = $ARGV[0] if $ARGV[0];
$job  = $ARGV[1] if $ARGV[1];

if ( (! $list) or (! $job) ) {
    die "Format: make_rhab_tars_14jun2021.pl  [input list of directories]  [batch job stem name] => print 1+ job scripts\n";
}

open my $LIST, '<', $list;

while (my $arch_dir = <$LIST>) {
    chomp $arch_dir;

    $i++;
    my $j = sprintf ("%02u", $i);

    my $outfile = "$job.$j.sh";
    $outfile    = safename($outfile);
    open my $OUTFILE, '>', $outfile;

    print $OUTFILE "#!/bin/bash\n";
    print $OUTFILE "#SBATCH --nodes=1\n";
    print $OUTFILE "#SBATCH --partition=RM-shared\n";
    print $OUTFILE "#SBATCH --time=48:00:00\n";
    print $OUTFILE "#SBATCH --ntasks-per-node=1\n";
    print $OUTFILE "#SBATCH --job-name=$outfile\n";
    print $OUTFILE "#SBATCH --mail-type=ALL\n";

    print $OUTFILE "cd $work_dir ;\n";
    
    print $OUTFILE "tar cvf $arch_dir.tar $arch_dir ;\n";
    print $OUTFILE "mv -i $arch_dir /ocean/projects/mcb190015p/schwarze/Z_backup/deprecated/rhabditella ;\n";
    print $OUTFILE "mv -i $arch_dir.tar /ocean/projects/mcb190015p/schwarze/Z_backup/tar_archives/rhabditella ;\n";

    close $OUTFILE;
}
close $LIST;

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

