#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

my $ofind   = q{};
my $es_list = q{};
my $label   = q{};

$ofind   = $ARGV[0] if $ARGV[0];
$es_list = $ARGV[1] if $ARGV[1];
$label   = $ARGV[2] if $ARGV[2];

my %es_homolog = ();
my @es_hits    = ();

if ( (! -e $ofind ) or (! -e $es_list ) or ( $label !~ /\A\S+\z/xms ) ) {
    die "Format: ofind2esp_06feb2023_v01.pl [OrthoFinder table] [ES gene list] [species label] > [ES homolog table]\n";
}

my $es_head = 'ES_' . $label;
my $header  = "Gene\t$es_head";

open my $ES_LIST, '<', $es_list;
while (my $esp = <$ES_LIST>) {
    chomp $esp;
    my $es_name = $esp . '(' . $label . ')';
    $es_homolog{$es_name} = 1;
}
close $ES_LIST;

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
        if ( exists $es_homolog{$name} ) {
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

