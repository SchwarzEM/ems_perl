#!/usr/bin/env perl

use strict;
use warnings;

my $list  = $ARGV[0];
my $table = $ARGV[1];

my %ok    = ();
my %valid = ();

open my $LIST, '<', $list or die "Can't open list of transposon PFAM motifs $list: $!";
while (my $input = <$LIST>) {
    chomp $input;
    if ( $input =~ /\A (PF\d+\.\d+) \z/xms ) { 
        my $motif   = $1;
        $ok{$motif} = 1;
    }
}
close $LIST or die "Can't close filehandle for list of transposon PFAM motifs $list: $!";

open my $TABLE, '<', $table or die "Can't open table of PFAM hits $table: $!";
while (my $input = <$TABLE>) {
    chomp $input;
    if ( $input =~ /\A (R[=]\d+) _\d+ \s .* \s (PF\d+\.\d+) \s/xms ) {
        my $element = $1;
        my $motif   = $2;
        if ( $ok{$motif}) { 
            $valid{$element} = 1;
        }
    }
}
close $TABLE or die "Can't close filehandle for table of PFAM hits $table: $!";

my @final_elements = sort keys %valid;
foreach my $final (@final_elements) {
    print "$final\n";
}


