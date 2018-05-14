#!/usr/bin/env perl

# update_wbg_ids.pl -- Erich Schwarz, 2/15/2010.
# Purpose: given a table with first 1 older and then 1 newer WBGene ID on each line, update names in a text stream.

use strict;
use warnings;
use Getopt::Long;

# Sample input (from TableMaker on WS210 gene IDs -- Gene ID followed by "merged to" Gene ID):
#
# "WBGene00016498"	"Caenorhabditis elegans"	"WBGene00003242"
#
# Another acceptable input would be:
# 
# WBGene00016498  WBGene00003242

my $wbid_table = q{};
my $old_id     = q{};
my $new_id     = q{};
my %old2new    = ();
my $help;

GetOptions ( 'table:s' => \$wbid_table,
             'help'    => \$help, );

if ( $help or (! $wbid_table ) ) { 
    die "update_wbg_ids.pl --table|-t [table of old-to-new WBGene IDs] <input files/stream>\n";
}

open my $WBID_TABLE, '<', $wbid_table or die "Can't open WormBase gene ID table file $wbid_table: $!"; 
while (my $input = <$WBID_TABLE>) { 
    chomp $input;

    # The word-boundary function will ignore quotations flanking names, e.g., those put by TableMaker.
    # But to be able to rely on end-of-line constraints, go ahead and strip out quotes anyhow.
    $input =~ s/\"//g;

    if ( $input =~ /\A \s* \b (WBGene\d+) \b .* \b (WBGene\d+) \b \s* \z /xms ) { 
        $old_id  = $1;
        $new_id  = $2;
        if ( exists $old2new{$old_id} ) { 
            die "Redundant mapping of old WBGene ID $old_id in: $input\n";
        }
        $old2new{$old_id} = $new_id;
    }

    else { 
        warn "Can't parse input line: $input\n";
    }
}
close $WBID_TABLE or die "Can't close filehandle to WormBase gene ID table file $wbid_table: $!";

while (my $input = <>) { 
    chomp $input;

    # Just going ahead and *maybe* replacing each old term -- blindly -- turns 
    #     out to be far more efficient than trying to scan the line for old terms.

    foreach my $obsolete_id (sort keys %old2new) { 
        my $correct_id = $old2new{$obsolete_id};
        $input =~ s/$obsolete_id/$correct_id/g;
    }

    print "$input\n";
}

