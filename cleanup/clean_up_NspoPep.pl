#!/usr/bin/perl

# basic_fasta_clean.pl, Erich Schwarz <emsch@its.caltech.edu>, 4/26/05.
# Purpose: clean up various irregularities in "FASTA" files.

while (<>) { 
    chomp ($input = $_);
    if ($input =~ /^>/) { 
        $input =~ m/^>\((.+?)\)(.+)/;
        $input = $1 . "  " . $2;
        while ($store) {
            $store =~ m/^(.{0,60})(.*)/;
            print "$1\n";
            $store = $2;
        }
        print ">"; print "$input\n";
        $input = "";
    }
    unless ($input =~ /^>/) {
        $input =~ s/[^a-zA-Z]//g;
        $store .= $input;
    }
}
while ($store) {
    $store =~ m/^(.{0,60})(.*)/;
    print "$1\n";
    $store = $2;
}
