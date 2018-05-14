#!/usr/bin/env perl

# assign_sORF_names.pl -- Erich Schwarz <emsch@its.caltech.edu> 11/5/2009.
# Purpose: provide application-friendly names to Shihan/Kou sORF CDSes.

use strict;
use warnings;
use Getopt::Long;

my $i      = 1;
my $suffix = q{};
my $help;

GetOptions ( 'suffix:s' => \$suffix,
             'help'     => \$help,   );

if ( ($help) or (! $suffix) ) { 
    die "Format: assign_sORF_names.pl",
        " --suffix|-s [suffix] <input stream/files>\n",
        ;
}
if ( $suffix and ( $suffix !~ / \A \w+ \z /xms ) ) { 
    die "Cannot use suffix $suffix!\n";
}

while (my $input = <>) { 
    chomp $input;
    if ( ( $input !~ /\A >/xms ) and ( $input =~ /\S/xms ) ) { 
        print "$input\n";
    }
    elsif ( $input =~ / \A > (\S .*) /xms ) { 
        my $headtext = $1;
        my $j = $i;
        $j = sprintf "%06d", $j;
        print '>sORF_', $suffix, '_', "$j   $headtext\n";
        $i++;
    }
    elsif ( $input =~ / \A > /xms ) { 
        die "Can't process: $input\n";
    }
}

