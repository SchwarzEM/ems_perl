#!/usr/bin/env perl

# altify_aug_preds.pl -- Erich Schwarz <ems394@cornell.edu>, 8/18/2020.
# Purpose: given AUGUSTUS gene names, stick 'alt' or some similar tag in front of their '.g1' etc.; useful when doing a merger of two AUGUSTUS predictions on the same genome.

use strict;
use warnings;
use autodie;
use Getopt::Long;

# Sample input line -- # start gene Raxei_Eur_2020.08.10.01.g1

my $infile   = q{};
my $pattern  = q{};
my $modifier = q{};
my $help;

GetOptions ( 'infile=s'   => \$infile,
             'pattern=s'  => \$pattern,
             'modifier=s' => \$modifier,
             'help'       => \$help,   );

$pattern ||= q{\S+?};

if ( $help or (! $infile) or (! $pattern) or (! $modifier) ) {
    die "Format: altify_aug_preds.pl\n",
        "    --infile|-i    [input file]\n",
        "    --pattern|-p   [name pattern recurring in header, e.g., ", q{"Raxei_Eur_2020\.08\.10.\d+"}, "; default is ", q{"\S+?"}, "]\n",
        "    --modifier|-m  [modifier such as 'alt' to stick in front of '.g1' etc. of original gene/transcript names]\n",
        "    --help|-h      [print this message]\n",
        ;
}

# Do special Perl quoting to make the pattern work better in a regex:
$pattern = qr/$pattern/;

open my $INPUT_FILE, '<', $infile;
while (my $input = <$INPUT_FILE>) { 
        chomp $input;
        $input =~ s/($pattern)(\.g\d+\w*)/$1.$modifier$2/xmsg;
        print "$input\n";
}
close $INPUT_FILE;


