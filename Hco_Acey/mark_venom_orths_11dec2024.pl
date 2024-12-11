#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $venom_list = q{};
my $orthologs  = q{};

$venom_list = $ARGV[0] if $ARGV[0];
$orthologs  = $ARGV[1] if $ARGV[1];

my %listed = ();

if ( (! $venom_list ) or (! $orthologs ) ) {
    die "Format: mark_venom_orths_11dec2024.pl [venom list with (taxon)s] [OrthoFinder table] > [OFind table with (taxon.venom) reassignments]\n";
}

open my $VENOM_LIST, '<', $venom_list;
while ( my $venom = <$VENOM_LIST> ) {
    chomp $venom;
    if ( $venom =~ /\A \S+ \z/xms ) {
        $listed{$venom} = 1;
    }
    else {
        die "From venom list $venom_list, cannot parse: $venom\n";
    }
}
close $VENOM_LIST;

open my $ORTHOLOGS, '<', $orthologs;
while ( my $input = <$ORTHOLOGS> ) {
    chomp $input;
    my $output = q{};
    while ( $input =~ /(\S+)(\s*)/xmsg ) {
        my $text  = $1;
        my $space = $2;
        if ( exists $listed{$text} ) {
            if ( $text =~ /\A \S+ \( \S+ \) \z/xms ) {
                $text =~ s/\)\z/.venom)/;
            }
            else {
                die "From orthologs file $orthologs, cannot reformat \"$text\" in: $input\n";
            }
        }
        $output = $output . $text . $space;
    }
    print "$output\n";
}
close $ORTHOLOGS;
