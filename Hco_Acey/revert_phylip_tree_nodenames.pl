#!/usr/bin/env perl

# revert_phylip_tree_nodenames.pl -- Erich Schwarz <ems394@cornell.edu>, 11/25/2013.
# Purpose: given a boostrapped PHYLIP tree with dorky names using reliable prefixes, the prefix, and a synonym table, recode the node names back from PHYLIP names to real ones.

use strict;
use warnings;

my $synonyms      = $ARGV[0];
my $phylip_prefix = $ARGV[1];
my $phylip_tree   = $ARGV[2];

my %phyl2real = ();

if (! ($synonyms and $phylip_prefix and $phylip_tree)) { 
    die "Format: revert_phylip_tree_nodenames.pl [synonym table] [phylip name prefix] [phylip-named tree]\n";
}

open my $SYNS, '<', $synonyms or die "Can't open synonym table $synonyms: $!";
while (my $input = <$SYNS>) { 
    chomp $input;

    # Sample input:
    # ASPphy_1        ASPR-s0042.g717/7-159
    # ASPphy_2        ASP-s0061.g3273/4-187

    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) { 
        my $phylip_name = $1;
        my $real_name   = $2;
        $phyl2real{$phylip_name} = $real_name;
    }
    else { 
        die "From synonym table $synonyms, can't parse input line: $input\n";
    }
}
close $SYNS or die "Can't close filehandle to synonym table $synonyms: $!";

open my $TREE, '<', $phylip_tree or die "Can't open PHYLIP tree $phylip_tree: $!";
while (my $input = <$TREE>) {
    chomp $input;
    
    # Sample input:
    # (((((((((((((((((ASPphy_64:100.0,ASPphy_65:100.0):100.0,((ASPphy_77:100.0,

    my $output      = q{};
    my $punctuation = q{};
    my $text        = q{};
    while ($input =~ / \A (.*?) ($phylip_prefix\d+) (.*) /xmsg ) {
        $punctuation = $1;
        $text        = $2;
        $input       = $3;
        $output .= $punctuation;
        if ( exists $phyl2real{$text} ) {
            $text = $phyl2real{$text};
        }
        $output .= $text;
    }
    $output .= $input;
    print "$output\n";
}
close $TREE or die "Can't close filehandle to PHYLIP tree $phylip_tree: $!";


