#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Cwd;
use File::Spec::Functions;  # catdir, catfile

my $working_dir = getcwd;

while (my $indir = <>) {
    chomp $indir;
    if ( $indir =~ /\A \S+ \/ ([^\s\/]+) \z/xms ) {
        my $target_subdir = $1;
        my $target_dir = catdir($working_dir, $target_subdir);
        print "mkdir -p $target_dir ;\n";
        print "rsync -av ",
              "$indir",
              '/* ',
              "$target_dir ;\n",
              ;
    }
    else {
        die "Cannot parse input: $indir\n";
    }
}

