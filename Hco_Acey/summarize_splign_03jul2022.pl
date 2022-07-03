#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Scalar::Util qw(looks_like_number);

my $total_nt   = 0;
my $total_id   = 0;
my $total_frac = 0;

while (my $input = <>) {
    chomp $input;
    my $fraction = 0;
    my $residues = 0;
    if ( $input =~ /\A (\S+) \t (\d+) \z/xms ) {
        $fraction = $1;
        $residues = $2;
        if ( looks_like_number($fraction) ) {
            $total_nt += $residues;
            $total_id += ($fraction * $residues);
        }
    }
}

if ( $total_nt >= 1 ) {
    $total_frac = ($total_id / $total_nt);
    $total_id   = int(($total_id + 0.5));
    $total_id   = commify($total_id);
    $total_nt   = commify($total_nt);
}

print "$total_frac\t$total_id\t$total_nt\n";


# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}

