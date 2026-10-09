#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);
use Cwd;

my $list = q{};
my $stem = q{};
my $dir  = getcwd;
my $i    = 0;

my @targets = ();

$list = $ARGV[0] if $ARGV[0];
$stem = $ARGV[1] if $ARGV[1];

if ( (! $list ) or (! $stem ) ) {
    die "Format: make_tar_w_checksum_08oct2026a.pl\n",
        "            [list of input dirs/files]\n",
        "            [batch job name stem]\n",
        "         => [N tar/checksum batch job files]\n",
        ;
}

open my $LIST, '<', $list;
while ( my $input = <$LIST> ) {
    chomp $input;
    if (! -e $input ) {
        die "From directory $dir, cannot find input dir/file: $input\n";
    }
    push @targets, $input;
}
close $LIST;

@targets = uniq(@targets);

foreach my $target (@targets) {
    $i++;
    my $j = sprintf "%02u", $i;
    my $outfile = 'job_'. "$stem.$j.sh";
    $outfile    = safename($outfile);

    my $tar_archive   = "$target.tar";
    my $tar_md5sum    = "$target.tar.md5";
    my $tar_sha256sum = "$target.tar.sha256";

    if ( ( -e $tar_archive ) or ( -e $tar_md5sum ) or ( -e $tar_sha256sum ) ) {
        die "Will not overwrite existing file: $tar_archive or $tar_md5sum or $tar_sha256sum\n";
    }

    open my $OUTFILE, '>', $outfile;

    print $OUTFILE '#!/bin/bash', "\n";
    print $OUTFILE '#SBATCH --nodes=1', "\n";
    print $OUTFILE '#SBATCH --partition=RM-shared', "\n";
    print $OUTFILE '#SBATCH --time=48:00:00', "\n";
    print $OUTFILE '#SBATCH --ntasks-per-node=4', "\n";
    print $OUTFILE "#SBATCH --job-name=$outfile\n";
    print $OUTFILE '#SBATCH --mail-type=ALL', "\n";
    print $OUTFILE "cd $dir ;\n";
    print $OUTFILE "tar cvf $tar_archive $target ;\n";
    print $OUTFILE "md5sum $tar_archive > $tar_md5sum ;\n";
    print $OUTFILE "sha256sum $tar_archive > $tar_sha256sum ;\n";

    close $OUTFILE;
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

