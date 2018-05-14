#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use File::Basename;

my @infiles = ();

while (my $input = <>) {
    chomp $input;
    push @infiles, $input;
}
        
print '#!/bin/bash', "\n";

foreach my $infile (@infiles) {
    my $dir = dirname($infile);
    my $file = basename($infile);
    my $check = "$infile.gz";
    if (! -r $check ) {
        die "Cannot read: $check\n";
    }
    print "cd $dir ;\n";
    print "ncftpput -u sra -p 'Qrjo6iJ4' ftp-private.ncbi.nih.gov / $file.gz ;\n";
}

