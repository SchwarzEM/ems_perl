#!/usr/bin/env perl

# rename_JPG.pl -- Erich Schwarz <emsch@its.caltech.edu>, 8/4/2008+.
# Purpose: rename *.JPGs in a directory to *.jpgs.
# 
# Note: this can be improved a LOT, e.g., with selectable directory and 
# suffixes.  But it's handy for simple work.

use strict;
use warnings;
use Cwd;

my $curr_dir = getcwd;

opendir my $DIR, $curr_dir 
    or die "Can't open current directory $curr_dir: $!";
my @files = readdir $DIR;

foreach my $file (@files) { 
    if (     ( $file !~ /\A \. /xms                 ) 
         and ( $file =~ / \A ( \S.* ) \.JPG \z /xms ) ) { 
        my $stem = $1;
        my $new_name = $stem . '.jpg';
        if (-e $new_name) {
            warn "Could not rename $file as $new_name,",   
                 " since file $new_name already exists!\n",
                 ;
        }
        if (! -e $new_name) { 
            rename $file, $new_name;
        }
    }
}

