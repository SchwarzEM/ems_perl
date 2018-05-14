#!/usr/bin/perl
# fasta2list.pl, Erich Schwarz <emsch@its.caltech.edu>, 4/24/05; v. simple.

while (<>) {
    $input = $_;
    if ($input =~ /^>(\S+)\s*/) {
        print "$1\n";
    }
}
