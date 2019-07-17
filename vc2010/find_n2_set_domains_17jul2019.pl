#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my %seen   = ();
my $header = q{};

my $set_table = q{};
my $n2_table  = q{};

$set_table = $ARGV[0] if $ARGV[0];
$n2_table  = $ARGV[1] if $ARGV[1];

if ( (! -r $set_table) or (! -r $n2_table ) ) { 
    die "Format: find_n2_set_domains_17jul2019.pl [SET domain gene table] [N2-to-VC2010 liftover table]\n";
}

open my $SET, '<', $set_table;
while (my $input = <$SET>) {
    chomp $input;
    if ( $input =~ /\A \S+ \t (chr\S+_pilon.g\d+) /xms ) {
        my $gene = $1;
        $seen{$gene} = $1;
    }
}
close $SET;

open my $N2, '<', $n2_table;
while (my $input = <$N2>) {
    chomp $input;
    if ( $input !~ /\S/xms ) { 
        die "Cannot parse empty text line from N2 table $n2_table\n";
    }
    elsif (! $header) {
        print "$input\n";
        $header = $input;
    }
    elsif ( $input =~ /\A \S+ \t [^\t]* \t [^\t]* \t (chr\S+_pilon.g\d+) /xms ) {
        my $gene = $1;
        if ( exists $seen{$gene} ) {
            print "$input\n";
        }
    }
}
close $N2;

