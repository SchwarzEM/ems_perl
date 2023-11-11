#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

my $ofind    = q{};
my $sub_list = q{};
my $subset   = q{};
my $species  = q{};

$ofind    = $ARGV[0] if $ARGV[0];
$sub_list = $ARGV[1] if $ARGV[1];
$subset   = $ARGV[2] if $ARGV[2];
$species  = $ARGV[3] if $ARGV[3];

my %sub_homolog = ();

if ( (! -e $ofind ) or (! -e $sub_list ) or ( $subset !~ /\A\S+\z/xms ) or ( $species !~ /\A\S+\z/xms ) ) {
    die "Format: ofind2subset_10nov2023.pl [OrthoFinder table] [subset gene list] [subset name] [species name] > [subset homolog table]\n";
}

my $sub_head = $subset . '_' . $species;
my $header   = "Gene\t$sub_head";

open my $SUB_LIST, '<', $sub_list;
while (my $subp = <$SUB_LIST>) {
    chomp $subp;
    my $sub_name = $subp . '(' . $species . ')';
    $sub_homolog{$sub_name} = 1;
}
close $SUB_LIST;

open my $OFIND, '<', $ofind;
while (my $input = <$OFIND>) {
    chomp $input;
    my $gene = q{};

    if ( $input =~ /\A (\S+) \b (.+) \z/xms ) {
        $gene  = $1;
        $input = $2;
    }
    else {
        die "From OrthoFinder $ofind, cannot parse: $input\n";
    }

    my @names = split '\s+', $input;
    my @homols = ();
    foreach my $name (@names) {
        if ( exists $sub_homolog{$name} ) {
            push @homols, $name;
        }
    }
    @homols = sort @homols;
    @homols = uniq @homols;
    my $homols_text = join '; ', @homols;
    if ( $homols_text =~ /\S/xms ) {
        print "$header\n" if $header;
        $header = q{};
        print "$gene\t$homols_text\n";
    }
}
close $OFIND;

