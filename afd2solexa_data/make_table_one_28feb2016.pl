#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

my @conds  = ();
my $header = "Conditional change\tGenes sig. up\tGenes sig. down";

my $data_ref;

while (my $input = <>) {
    chomp $input;
    if ( $input =~ /\A (.+ vs\. .+) [ ] (\S+)reg\. \t (\d+) \z/xms ) { 
        my $cond  = $1;
        my $dir   = $2;
        my $count = $3;
        push @conds, $cond;
        if ( ( $dir ne 'up' ) and ( $dir ne 'down' ) ) {
            die "Cannot parse direction of change in: $input\n";
        }
        $data_ref->{'cond'}->{$cond}->{'dir'}->{$dir} = $count;
    }
    elsif (     ( $input !~ /\A \s* \z/xms ) 
            and ( $input !~ /\A Conditional [ ] change \t Genes [ ] w\/ [ ] sig\. [ ] change \z/xms ) ) { 
        die "Cannot parse input line: $input\n";
    }
}

@conds = uniq(@conds);

foreach my $cond (@conds) {
    my $up_count   = $data_ref->{'cond'}->{$cond}->{'dir'}->{'up'};
    my $down_count = $data_ref->{'cond'}->{$cond}->{'dir'}->{'down'};
    $up_count      = commify($up_count);
    $down_count    = commify($down_count);
    print "$header\n" if $header;
    $header = q{};
    print "$cond\t$up_count\t$down_count\n";
}

# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}


