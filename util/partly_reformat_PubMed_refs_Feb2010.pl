#!/usr/bin/env perl

# partly_reformat_PubMed_refs.pl -- Erich Schwarz <emsch@its.caltech.edu>, 12/7/2009.
# Purpose: convert (e.g) "Lin JJ, Lobo NF, Lopez JR, Malek JA" to "Lin, J.J, Lobo, N.F., Lopez, J.R., Malek, J.A.".

use strict;
use warnings;

my $input_text     = q{};
my $remaining_text = q{};
my $output_text    = q{};

while (my $input = <>) {
    chomp $input;
    $input_text .=  q{ } . $input;
}

while ( $input_text =~ /\A ([^,]+) \s ([A-Z]+) \, (.+) /xms ) {
    my ($surname, $initials);
    $surname        = $1;
    $initials       = $2;
    $remaining_text = $3;
    $output_text .= "$surname, ";
    while ( $initials =~ / ([A-Z]{1}) (.*) /xms ) {
        my $front_text = $1;
        $initials = $2;
        $output_text .= "$front_text.";
        $initials = $2;
    }
    $output_text .= q{,};
    $input_text = $remaining_text;
}
$output_text .= $remaining_text;

$output_text =~ s/\A\s+//;
$output_text =~ s/\s+\z//;

# Export as one line, but spaced to make cut-and-paste easy.
print "\n";
print "$output_text\n";
print "\n";

