#!/usr/bin/env perl

# make_md5sum_file.pl -- Erich Schwarz <emsch@its.caltech.edu>, 10/5/2010.
# Purpose: given files in a directory, print out their names and md5sum values.

use strict;
use warnings;
use Cwd;
use File::Basename;

my $md5sum_text = q{};
my $curr_dir = getcwd;

opendir my $DIR, $curr_dir 
    or die "Can't open current directory $curr_dir: $!";
my @files = readdir $DIR;

foreach my $file (@files) { 
    if (-f $file ) { 
        $md5sum_text = `md5sum $file`;
        chomp $md5sum_text;
        print "$md5sum_text\n";
    }
}

