#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use autodie;

my $data_ref;

my $wbgenes  = q{};
my $families = q{};
my $some;
my $reject;
my $help;

GetOptions ( 'wbgenes=s'  => \$wbgenes,
             'families=s' => \$families,
             'some'       => \$some,
             'reject'     => \$reject,
             'help'       => \$help,
);

if ( $help or (! $wbgenes) or (! $families) ) {
    print "Format: select_fams_w_wbgene_member_set.pl\n",
        "    --wbgenes|-w   [List of WBGene IDs, one ID per line, which family members will be required to have, by default exclusively]\n",
        "    --families|-f  [List of gene families, one family per line, which will be checked for having WBGene IDs from the list (solely, by default; options change this)]\n",
        "    --some|-s      [Optionally, print gene families which have *any* members from the wbgenes list, even if they also have non-members]\n",
        "    --reject|-r    [Optionally, reject gene families with members from the wbgenes list, by default if their members solely come from the list;\n",
        "                       if used in conjunction with the --some argument, this rejects gene families with *any* member from the wbgenes list]\n",
        "    --help|-h      [print this message]\n",
        ;
    exit;
}

open my $WBGENES, '<', $wbgenes;
while ( my $input = <$WBGENES>) { 
    chomp $input;
    if ( $input =~ / WBGene\d+ .* WBGene\d+ /xms ) {
        die "In WBGene list file $wbgenes, two or more IDs per line, in: $input\n";
    }
    elsif ( $input =~ / (WBGene\d+) /xms ) { 
        my $wbgene_id = $1;
        $data_ref->{'listed_wbgene_id'}->{$wbgene_id} = 1;
    }
    else {
        die "From WBGene list file $wbgenes, can't parse: $input\n";
    }
}
close $WBGENES;

open my $FAMILIES, '<', $families;
while ( my $input = <$FAMILIES>) {
    chomp $input;

    # For each line, initialize these array variables as zero.
    my @wb_gene_ids_listed     = ();
    my @wb_gene_ids_not_listed = ();

    while ( $input =~ / (WBGene\d+) /xmsg ) {
        my $wbgene_id = $1;
        if ( exists $data_ref->{'listed_wbgene_id'}->{$wbgene_id} ) {
            push @wb_gene_ids_listed, $wbgene_id;
        }
        else {
            push @wb_gene_ids_not_listed, $wbgene_id;
        }
    }

    # All-or-nothing: we either print lines with exclusively listed membership, or lines with exclusively non-listed membership.
    if (! $some) {
        if ( (! $reject) and @wb_gene_ids_listed and (! @wb_gene_ids_not_listed) ) {
            print "$input\n";
        }
        elsif ( $reject and (! @wb_gene_ids_listed) and @wb_gene_ids_not_listed ) {
            print "$input\n";
        }
    }

    # Conditional behavior: we either print lines with *some* listed membership, or lines with *some* non-listed membership.
    elsif ($some) {
        if ( (! $reject) and @wb_gene_ids_listed ) {
            print "$input\n";
        }
        elsif ( $reject and @wb_gene_ids_not_listed ) {
            print "$input\n";
        }
    }
}
close $FAMILIES;


