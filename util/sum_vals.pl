#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Scalar::Util qw(looks_like_number);

my $sum = 0;

while (my $input = <>) {
    chomp $input;

    $input =~ s/\A\s+//;
    $input =~ s/\s+\z//;

    if (! looks_like_number($input) ) {
        die "Input text not a number: $input\n";
    }
    $sum = $sum + $input;
}

$sum = sprintf("%.2f", $sum);
$sum = commify($sum);

print "$sum\n";

# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}

