#!/usr/bin/perl

# clean_TIGR_rice.pl, Erich Schwarz <emsch@its.caltech.edu>, 7/17/05.
# Purpose: clean up various irregularities in TIGR rice protein FASTA from
#     ftp://ftp.tigr.org/pub/data/Eukaryotic_Projects/o_sativa/annotation_dbs/pseudomolecules/version_3.0/all_chrs/all.seq

while (<>) { 
    chomp ($input = $_);
    if ($input =~ /^>/) { 
        while ($store) {
            $store =~ m/^(.{0,60})(.*)/;
            print "$1\n";
            $store = $2;
        }
        $input =~ s/\|/  /g;
        print "$input\n";
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
