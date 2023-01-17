#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $listfile = q{};
my $genbank  = q{};

$listfile = $ARGV[0] if $ARGV[0];
$genbank  = $ARGV[1] if $ARGV[1];

my %ok_gb = ();
my $print = 0;

open my $LISTFILE, '<', $listfile;
while (my $input = <$LISTFILE>) {
    chomp $input;
    if ( $input =~ /\A \S+ \z/xms ) {
        $ok_gb{$input} = 1;
    }
    else {
        die "From list file $listfile, cannot parse: $input\n";
    }
}
close $LISTFILE;

open my $GENBANK, '<', $genbank;
while (my $input = <$GENBANK>) {
    chomp $input;
    if ( $input =~ /\A LOCUS \s+ (\S+) /xms ) {
        my $locus = $1;
        if ( exists $ok_gb{$locus} ) {
            $print = 1;
        }
        else {
            $print = 0;
        }
    }
    elsif ( $input =~ / LOCUS /xms ) {
        die "From GenBank file $genbank, cannot parse: $input\n";
    }
    if ($print) {
        print "$input\n";
    }
}
close $GENBANK;


