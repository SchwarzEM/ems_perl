#!/usr/bin/perl

# single_line_cdhit_clustr_archival.pl -- Erich Schwarz <emsch@its.caltech.edu>, sometime in the Pleistocene.
# Purpose: convert *.clstr outputs of CD-HIT into single-line summaries; this version is ARCHIVAL, and strongly discouraged.  Use the modern one!

$read_line = 0;
while (<>) {
    chomp ($input = $_);
    if ($input =~ /^>/) { 
        if ($read_line) { print "\n"; }
    }
    elsif ($input =~ /[^>]*>(\S+)(.*)/) {
        $read_line = 1;
        while ($input =~ /[^>]*>(\S+)(.*)/) {
            print "$1\t";
            $input = $2;
        }
    }
}
print "\n";
