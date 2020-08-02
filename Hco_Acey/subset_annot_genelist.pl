#!/usr/bin/env perl

# subset_annot_genelist.pl -- Erich Schwarz <ems394@cornell.edu>, 8/2/2020.
# Purpose: list A has good gene order; list B has wanted genes in bad order; read both, then output B's genes in A's order.  Requires that A >= B.

use strict;
use warnings;
use autodie;

my $well_ordered_list   = q{};
my $wanted_members_list = q{};

$well_ordered_list   = $ARGV[0] if $ARGV[0];
$wanted_members_list = $ARGV[1]	if $ARGV[1];

if ( (! $well_ordered_list) or (! $wanted_members_list) ) {
    die "Format: subset_annot_genelist.pl [well-ordered list] [wanted-members list] > [wanted members, well-ordered]\n";
}

my %selected = ();

open my $WANT, '<', $wanted_members_list;
while (my $input = <$WANT>) {
    chomp $input;
    $input =~ s/\A\s+//;
    $input =~ s/\s+\z//;
    if ( $input =~ /\S/xms ) { 
        $selected{$input} = 1;
    }
}
close $WANT;

open my	$ORDER, '<', $well_ordered_list;
while (my $input = <$ORDER>) {
    chomp $input;
    if ( exists $selected{$input} ) {
        print "$input\n";
    }
}
close $ORDER;

