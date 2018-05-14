#!/usr/bin/perl

# rename_files_07jul2006.pl -- Erich Schwarz <emsch@its.caltech.edu>, 7/7/2006.
# Purpose: do simple renaming of non-Unixish files.

use strict;
use warnings;

chdir $ARGV[0];

my $oldfile = "";
my $newfile = "";
my @files = ();
opendir (MYDIR, ".");
@files = readdir MYDIR;
foreach $oldfile (@files) { 
    if (-f $oldfile) { 
           $newfile = $oldfile;
           $newfile =~ s/[^a-zA-Z0-9\.]/_/g;
           $newfile =~ s/_+/_/g;
           $newfile =~ s/\._/_/g;
           if (-s $newfile) { print "Can't rename:\t$oldfile\t$newfile!\n"; }
           unless (-s $newfile) { 
               rename $oldfile, $newfile;
               print "Renamed:\t$oldfile\t$newfile.\n";
           }
           $newfile = "";
        }
}
