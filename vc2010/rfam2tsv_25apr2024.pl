#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my @data_col_names = (
    "idx",
    "target name",
    "accession",
    "query name",
    "accession",
    "clan name",
    "mdl",
    "mdl from",
    "mdl to",
    "seq from",
    "seq to",
    "strand",
    "trunc",
    "pass",
    "gc",
    "bias",
    "score",
    "E-value",
    "inc",
    "olp",
    "anyidx",
    "afrct1",
    "afrct2",
    "winidx",
    "wfrct1",
    "wfrct2",
    "mdl len",
    "seq len",
    "description of target",
);

my $header = join "\t", @data_col_names;

while ( my $input = <> ) {
    chomp $input;
    if ( $input !~ /\A[#]/xms ) {
        $input =~ s/\s+\z//;

        # This precise typography is required to capture $front properly: '((?:\S+ \s+){28})'
        if ( $input =~ /\A ((?:\S+ \s+){28}) (.+) \z/xms ) {
            my $front = $1;
            my $desc  = $2;

            $front =~ s/\s+\z//;
            $desc  =~ s/\s+\z//;

            my @vals = split /\s+/, $front;
            my $text = join "\t", @vals;
            $text    = "$text\t$desc";

            print "$header\n" if $header;
            $header = q{};

            print "$text\n";
        }
        else {
            die "Cannot parse: $input\n";
        }
    }
}
