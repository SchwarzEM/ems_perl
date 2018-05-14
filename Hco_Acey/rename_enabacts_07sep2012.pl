#!/usr/bin/env perl

use strict;
use warnings;

my %seen = ();

while (my $input = <>) {
    chomp $input;
    if ( $input !~ /\A \> /xms ) { 
        print "$input\n";
    }
    elsif ( $input =~ /\A \> (\S+\s+ (\S+) \s.+) \z/xms ) { 
        my $orig_header = $1;
        my $accession   = $2;
        if ( exists $seen{$accession} ) { 
            $accession = safename_id($accession,\%seen);
            if ( exists $seen{$accession} ) {
                die "Subroutine safename_id failed to suffix $accession into original name\n";
            }
            warn "Forced to rename sequence to $accession\n";
        }
        $seen{$accession} = 1;
        print '>BACT_', "$accession $orig_header\n", ;
    }
    else { 
        die "Can't parse input line: $input\n";
    }
}

sub safename_id { 
    my $_acc_orig     = $_[0];
    my $_seen_hashref = $_[1];
    my $_acc          = $_acc_orig;
    my $i = 1;
    while ( exists $_seen_hashref->{$_acc} ) {
        $i++;
        $_acc = $_acc_orig . ".$i";
    }
    return $_acc;
}

