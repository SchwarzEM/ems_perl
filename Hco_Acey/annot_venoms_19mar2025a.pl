#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $venom_list  = q{};
my $gene_annots = q{};

$venom_list  = $ARGV[0] if $ARGV[0];
$gene_annots = $ARGV[1] if $ARGV[1];

if ( (! $venom_list ) or (! $gene_annots ) ) {
    die "Format: annot_venoms_19mar2025a.pl [venom list] [gene annots] > [venom annots]\n";
}

my %venoms = ();

open my $VENOM_LIST, '<', $venom_list;
while ( my $input = <$VENOM_LIST> ) {
    chomp $input;
    if ( $input =~ /\A \S+ \z/xms ) {
        $venoms{$input} = 1;
    }
    else {
        die "From venom list $venom_list, cannot parse: $input\n";
    }
}
close $VENOM_LIST;

open my $GENE_ANNOTS, '<', $gene_annots;
while ( my $input = <$GENE_ANNOTS> ) {
    chomp $input; 
    if ( $input =~ /\A (\S+) \t ([^\t]+) \z/xms ) {
        my $gene    = $1;
        my $id_text = $2;
        my @ids     = split '; ' , $id_text;
        my @matches = ();
        if ( $gene eq 'Gene' ) {
            print "Gene\tVenom_match\n";
        }
        else {
            foreach my $id (@ids) {
                if ( exists $venoms{$id} ) {
                    my $match = "Venom [$id]";
                    push @matches, $match;
                 }    
            }
            if ( @matches ) {
                my $match_text = join '; ', @matches;
                print "$gene\t$match_text\n";
                @matches = ();
            }
       }
    }
    else {
        die "From gene annots file $gene_annots, cannot parse; $input\n";
    }
}
close $GENE_ANNOTS;
