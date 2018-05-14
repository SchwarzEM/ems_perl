#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $good_list = $ARGV[0];
my $mauve_dat = $ARGV[1];

my %good_contigs = ();

my %ori_code = (
    'forward'    => '+',
    'complement' => '-',
);

my $read_lines = 0;


open my $GOOD, '<', $good_list;
while (my $input = <$GOOD>) {
    chomp $input;

    # sample desired input line:
    # >nigoni_2015.12.01_059  352,027 nt      1K      5K

    if ( $input =~ /\A [>] (\S+) .* \s 5K \s* \z/xms ) { 
        my $contig = $1;
        $good_contigs{$contig} = 1;
    }
}
close $GOOD;

open my $MAUVE, '<', $mauve_dat;
while (my $input = <$MAUVE> ) {
    chomp $input;

    # sample desired input line:
    # contig  nigoni_2015.12.01_014   chromosome      forward 1       3241740

    if ( $input =~ /\A Ordered \s+ Contigs \s* \z/xms ) { 
        $read_lines = 1;
    }
    elsif ( $input =~ /\A Contigs \s+ with \s+ conflicting \s+ ordering \s+ information \s* \z/xms ) { 
        $read_lines = 0;
    }
    elsif ( $read_lines and ( $input =~ /\A contig \s+ (\S+) \s+ chromosome \s+ (forward|complement) \s+ \d+ \s+ \d+ \s* \z/xms ) ) { 
        my $contig = $1;
        my $orient = $2;
        $orient = $ori_code{$orient};
        if ( exists $good_contigs{$contig} ) {
            print "$contig\t$orient\n";
        }
    }
}

