#!/usr/bin/env perl

# sieve_WBPmat_ace.pl -- Erich Schwarz <emsch@its.caltech.edu>, 9/25/2008.
# Purpose: given a .ace file with PFM and/or PWM Position_Matrix objects, echo-print only PFMs.

use strict;
use warnings;

my $type   = shift @ARGV;
my $object = q{};
my %valid  = ();

my $object2info_ref;

if ( ($type ne 'Frequency' ) and ($type ne 'Weight' ) ) { 
    die "Format: ./sieve_PFM_ace.pl Frequency|Weight [*.ace]\n";
}

while (my $input = <>) { 
    if ( $input =~ /\A \s* \z/xms ) { 
        $object = q{};
    }
    # Order of next two matters; reverse order double-print starting line!
    if ( ($object) and ( $input =~ /\S/xms ) ) {
        if ( $input =~ /\A Type \s+ $type /xms ) {
            $valid{$object} = 1;
        }
        push @{ $object2info_ref->{$object} }, $input;
    }
    if ( (! $object) and ( $input =~ /\S/xms ) ) { 
        if ( $input =~ /\A Position_Matrix \s+ : \s+ \" (WBPmat\d+) \"  /xms ) { 
            $object = $1;
            push @{ $object2info_ref->{$object} }, $input;
        }
    }
}

foreach my $ok_object (sort keys %valid) { 
    print "\n";
    foreach my $line ( @{ $object2info_ref->{$ok_object} } ) { 
        print $line;
    }
    print "\n";
}
