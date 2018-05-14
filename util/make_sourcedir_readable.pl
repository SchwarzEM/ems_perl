#!/usr/bin/perl

# make_sourcedir_readable.pl
# Erich Schwarz <emsch@its.caltech.edu>, 9/12/04
# Purpose: make excessively-secure (root umask 077) source trees (e.g. Perl lib) readable.

use File::Find;

@ARGV = qw(.) unless @ARGV;

find (\&simple, @ARGV);

sub simple { 
    if (-d $_) {
        chmod 0755, $_; 
        print "$File::Find::name\n"; 
    }
    elsif (-f $_) {
        chmod 0644, $_;
        print "$File::Find::name\n";
    }
}
