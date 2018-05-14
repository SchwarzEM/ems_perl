#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $ofinder = q{};
my $sp1     = q{};
my $sp2     = q{};
my $prefix  = q{};

$ofinder = $ARGV[0] if $ARGV[0];
$sp1     = $ARGV[1] if $ARGV[1];
$sp2     = $ARGV[2] if $ARGV[2];
$prefix  = $ARGV[3] if $ARGV[3];

if (    (! -r $ofinder ) 
     or ( $sp1 !~ /\A \S+ \z/xms ) 
     or ( $sp2 !~ /\A \S+ \z/xms ) 
     or ( $prefix !~ /\A \S+ \z/xms ) 
   ) {
    die "Format: ncbify_strict_ofinder_pair_12oct2017.pl",
        " [OrthoFinder file] [species 1, 1:1 query] [species 2, 1:1 target] [prefix]",
        " > [output for NCBI/GenBank gene synonyms]\n",
        ;
}

open my $OFINDER, '<', $ofinder;
while (my $input = <$OFINDER>) {
     chomp $input;
     if ( $input =~ /\A [^\t]+ \t ([^\t]+) \z/xms ) { 
         my $orthodata = $1;

         my @taxa = ($sp1, $sp2);
         foreach my $taxon (@taxa) {
             if ( $orthodata =~ / \($taxon\) .* \($taxon\) /xms ) { 
                 die "In OrthoFinder file $ofinder, observe 2+ instances of $taxon in: $input\n";
             }
             if ( $orthodata !~ / \b \S+ \($taxon\) /xms ) {
                 die "In OrthoFinder file $ofinder, fail to observe any usable instance of $taxon in: $input\n";
             }
         }

         my $sp1_gene = q{};
         if ( $orthodata =~ / \b (\S+) \($sp1\) /xms ) {
             $sp1_gene = $1;
         }
         else {
             die "No! it cannot be! -- $input\n";
         }

         my $sp2_gene = q{};
         if ( $orthodata =~ / \b (\S+) \($sp2\) /xms ) {
             $sp2_gene = $1;
         }
         else {
             die "So near, and yet so far. -- $input\n";
         }
         
         # trim down extra namecruft
         if ( $sp1_gene =~ /\A \S+ \| (\S+?) \z/xms ) {
             $sp1_gene = $1;
         }
         if ( $sp2_gene =~ /\A \S+ \| (\S+?) \z/xms ) {
             $sp2_gene = $1;
         }

         # finally:
         print "$sp1_gene\tgene_synonym\t$prefix$sp2_gene\n";
     }
     else { 
         die "In OrthoFinder file $ofinder, cannot parse: $input\n";
     }
}
close $OFINDER;
