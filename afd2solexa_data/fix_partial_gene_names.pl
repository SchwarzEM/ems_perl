#!/usr/bin/env perl

# fix_partial_gene_names.pl -- Erich Schwarz <emsch@its.caltech.edu>, 5/26/2011.
# Purpose: take a name-crippled file and turn it back into a solid one.

use strict;
use warnings;

my $name_table = $ARGV[0];
my $broken_file = $ARGV[1];

my ($wbgene, $cds, $cgc, $name, $rest_of_line);
my %name2gene = ();

open my $NAMES, '<', $name_table or die "Can't open name table $name_table: $!";
while (my $input = <$NAMES>) {
    chomp $input; 
    if ( $input =~ / \A (WBGene\d+) \t ([^\t]+) \t ([^\t]+) /xms ) { 
        $wbgene = $1;
        $cds    = $2;
        $cgc    = $3;
        $wbgene = $wbgene . '|' . $cds . '|' . $cgc;
        $name2gene{$cds} = $wbgene;
        $name2gene{$cgc} = $wbgene;
    }
    elsif ( $input =~ / \A (WBGene\d+) \t ([^\t]+) \t /xms ) {
        $wbgene = $1;
        $cds    = $2;
        $wbgene = $wbgene . '|' . $cds ;
        $name2gene{$cds} = $wbgene;
    }
    else { 
        die "Can't parse input line from name table $name_table: $input!\n";
    }
}
close $NAMES or die "Can't close filehandle to name table $name_table: $!";

open my $BROKEN, '<', $broken_file or die "Can't open name-broken file $broken_file: $!";
while (my $input = <$BROKEN>) {
    chomp $input;
    if ( $input =~ / \A (\S+) (.*) \z /xms ) { 
        $name         = $1;
        $rest_of_line = $2;
        if ( exists $name2gene{$name} ) {
            $input = $name2gene{$name} . $rest_of_line;
        }
    }
    print "$input\n";
}
close $BROKEN or die "Can't close filehandle to name-broken file $broken_file: $!";

