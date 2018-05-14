#!/usr/bin/perl

# strip_lowercase_dna.pl, Erich Schwarz <emsch@its.caltech.edu>, 7/20/05.
# Purpose: strip out lowercase nucleotides, output result as clean FASTA.

while (<>) {
    chomp ($input = $_);
    if ($input =~ /^>/) {
        while ($store) {
            $store =~ m/^(.{0,60})(.*)/;
            print "$1\n";
            $store = $2;
        }
        print "$input\n";
    }
    unless ($input =~ /^>/) {
        $input =~ s/[^a-zA-Z]//g;
        $input =~ s/[a-z]//g;
        $store .= $input;
    }
}
while ($store) {
    $store =~ m/^(.{0,60})(.*)/;
    print "$1\n";
    $store = $2;
}
