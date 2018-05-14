#!/usr/bin/env perl

# part2wbfullgenes.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/16/2010.
# Purpose: given an incoming CDS or WBGeneID list stream and a defined wormpep, emit full ID stream.

use strict;
use warnings;
use Getopt::Long;

my $names  = q{};

my $help;
my $pretty;

my %part2fullname = ();
my %possible_part = ();

GetOptions ( 'names=s', => \$names,
             'help'     => \$help,  
             'pretty'   => \$pretty, );

# Typical wormpep190 lines:
# >4R79.1b        CE39659 WBGene00003525  locus:nas-6 ..
# >AC7.3  CE07653 WBGene00014997

if ( $help or (! $names) ) {   
    die "Format: part2wbfullgenes.pl  --names|-n [WBGene/cds/cgc table, or wormpep]",
        "  (optional: --pretty|-p)",
        "  [files or STDIN]\n",
        ;
}

open my $NAMES, '<', $names 
    or die "Can't open WBGene/cds/cgc table or wormpep ($names): $!";
while (my $input = <$NAMES>) { 
    chomp $input;
    # Read names from wormpep:
    if ( $input =~ /\A > (\S+) \s+ CE\d+ \s+ (WBGene\d+) \s+ (\S+) /xms ) { 
        my $cds         = $1;
        my $gene        = $2;
        my $maybe_locus = $3;
        my $locus       = q{};
        $cds =~ s/[a-z]\z//;
        my @id_tags = ($gene, $cds);
        if ( $maybe_locus =~ /\A locus : (\S+) \z/xms ) { 
            $locus = $1;
            push @id_tags, $locus;
        }
        my $full_id = join '|', @id_tags;
        $part2fullname{$cds}  = $full_id;
        $part2fullname{$gene} = $full_id;
    }

    # Or, read names from plain tab-delimited WBGene/CDS/CGC table:
    elsif ( $input =~ /\A (WBGene\d+) \s+ (\S.+\S) \s* \z/xms ) { 
        my $gene        = $1;
        my $other       = $2;
        my $cds         = q{};
        my $locus       = q{};
        my @id_tags = ($gene);
        if ( $other =~ /\A (\S+) \s+ (\S+) \z/xms ) {
            $cds   = $1;
            $locus = $2;
            push @id_tags, $cds, $locus;
        }
        elsif ( $other =~ /\A (\S+) \z/xms ) { 
            $cds   = $1;
            push @id_tags, $cds;
        }
        else { 
            die "Can't parse \"$other\" in \"$input\"!\n";
        }
        my $full_id = join '|', @id_tags;
        $part2fullname{$cds}  = $full_id;
        $part2fullname{$gene} = $full_id;
    }
}
close $NAMES 
    or die "Can't close filehandle to WBGene/cds/cgc table or wormpep ($names): $!";

while (my $input = <>) { 
    chomp $input;
    my $output = $input;
    # Use non-greedy undefined residues, so that lines with *no* "\s" following the ID can be parsed:
    if ( $input =~ /\A (\S+) (.*?) \z /xms ) { 
        my $partial_id   = $1;
        my $rest_of_line = $2;
        if ( $part2fullname{$partial_id} ) { 
            my $printed_name = $part2fullname{$partial_id};
            if ($pretty) { 
                $printed_name = sprintf "%-35s", $printed_name;
            }
            $output = $printed_name . $rest_of_line;
        }
        if (! $part2fullname{$partial_id} ) { 
            warn "Can't parse: $input\n";
        }
    }
    print "$output\n";
}

