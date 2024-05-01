#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my %map = ();

my $ids2terms = q{};
my $go_annot  = q{};

$ids2terms = $ARGV[0] if $ARGV[0];
$go_annot  = $ARGV[1] if $ARGV[1];

if ( (! $ids2terms) or (! $go_annot) ) {
    die "Format: expand_goterms_30apr2024.pl [ids2terms] [go_annot] > [term-expanded go_annot] ;\n";
}

open my $IDS2TERMS, '<', $ids2terms;
while ( my $input = <$IDS2TERMS> ) {
    chomp $input;
    if ( $input =~ /\A (GO:\d+) \t (.*) \z/xms ) {
        my $id   = $1;
        my $term = $2;
        $term =~ s/\A\s+//;
        $term =~ s/\s+\z//;
        if ( $term !~ /\S/xms ) {
            die "In $ids2terms, cannot map GO ID $id to useful term (\"$term\")\n";
        }
        if ( exists $map{$id} ) {
            die "In $ids2terms, redundant mapping of GO ID $id to both $map{$id} and $term\n";
        }
        $map{$id} = $term;
    }
}
close $IDS2TERMS;

open my $GO_ANNOT, '<', $go_annot;
while ( my $input = <$GO_ANNOT> ) {
    chomp $input;
    $input =~ s/(GO:\d+)/$map{$1} ($1)/g;
    print "$input\n";
}
close $GO_ANNOT;

