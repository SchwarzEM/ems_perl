#!/usr/bin/perl
# fastaheads.pl, Erich Schwarz <emsch@its.caltech.edu>, 1/18/06; v. simple.

while (<>) {
    chomp($input = $_);
    if ($input =~ /^>(\S.*)/) {
        print "$1\n";
    }
}
