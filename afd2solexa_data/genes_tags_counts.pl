#!/usr/bin/perl

# genes_tags_counts.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/28/2008.
# Purpose: get useful gene/tag counts from Solexa splice/uniqs.

use strict;
use warnings;

# Enforce superficially correct inputs.
my $splice_file = $ARGV[0];
my $uniqs_file  = $ARGV[1];

unless ( $#ARGV == 1 ) {
    die 'Format: ./genes_tags_counts.pl  [splice file]  [uniqs file]', "\n";
}

# These hashes both record a sighting and count tags.
my %genes_spliced = ();
my %genes_tagged  = ();
my %ambivs_tagged = ();

# Input lines of uniqs vs. splices are ~same, so ~same parsing:
sub parse_input_line { 
    my ($_input, $_splice);
    ($_input, $_splice) = @_;

    # Typical input lines:
    # WBGene00019904|WBGene00020894   42956 -- never in splices, though.
    # WBGene00020894  13716
    # WBGene00007790  1

    if ( $_input =~ / \A ( WBGene\d+ (?: \| WBGene\d+ )* )
                            \s+
                            ( \d+ )  
                         \z  /xms ) {
        my $gNtag = $1;
        my $count = $2;

        # Splices are stronger evidence, so track them specifically:
        if ($_splice eq "splice") { 
            $genes_spliced{$gNtag} += $count;
        }

        # Genes can be ambiguous, because tags can touch 2 exons:
        if ( $gNtag =~ / \| /xms ) {
            $ambivs_tagged{$gNtag} += $count;
        }
        if ( $gNtag !~ / \| /xms ) {
            $genes_tagged{$gNtag}  += $count;
        }
    }
}

# Import data from the splices and uniqs files.

open my $SPLICES, "<", "$splice_file" 
    or die "Can't open splice file $splice_file\n";

while (my $input = <$SPLICES>) {
    chomp $input;
    parse_input_line($input, "splice");
}

close $SPLICES;

open my $UNIQS, "<", "$uniqs_file"
    or die "Can't open uniqs. file $uniqs_file\n";

while (my $input = <$UNIQS>) { 
    chomp $input;
    parse_input_line($input, q{});
}

# Export a summary of uniqs/splices data.

# A note on sorting genes in the summary:
# 
# Bare "sort" (or 'cmp') is used to give clean ASCIIbetical order.
# $b <=> $a is then used to give genes in descending count order.
# And int($b <=> $a) comparisons are used to avoid floating-point nonsense.
# So a cohort of equal-counted genes will be listed reliably ASCIIbetically.

my @Ntagged_genes 
    = sort { int($genes_tagged{$b}) <=> int($genes_tagged{$a}) } 
      sort 
      grep { int($genes_tagged{$_}) != 1 } 
      keys %genes_tagged;
my $tagged_count = @Ntagged_genes;

my @Ntagged_ambivs 
    = sort { int($ambivs_tagged{$b}) <=> int($ambivs_tagged{$a}) } 
      sort 
      grep { int($ambivs_tagged{$_}) != 1 } 
      keys %ambivs_tagged;
my $ambiv_count = @Ntagged_ambivs;

my @unitagged_splices 
    = sort
      grep { int($genes_tagged{$_})  == 1 }
      grep { int($genes_spliced{$_}) == 1 }
      keys %genes_spliced;
my $uni_splic_count = @unitagged_splices;

my @unitagged_genes  
    = sort 
      grep {! $genes_spliced{$_} }
      grep { int($genes_tagged{$_}) == 1 } 
      keys %genes_tagged;
my $unitagged_count = @unitagged_genes;

my @unitagged_ambivs 
    = sort 
      grep { int($ambivs_tagged{$_}) == 1 } 
      keys %ambivs_tagged;
my $uni_ambiv_count = @unitagged_ambivs;

print "Genes with >= 2 tags:                $tagged_count\n";
print "Ambiv. id. genes with >= 2 tags:     $ambiv_count\n";
print "Genes with 1 splice tag:             $uni_splic_count\n";
print "Genes with 1 non-splice tag:         $unitagged_count\n";
print "Ambiv. genes w/ 1 non-splice tag:    $uni_ambiv_count\n";

print "\n";

# There are a number of ways one might list the genes.
# One reasonable way, though is to have a straight list
#    of both unambig. and ambig. genes with >= 2 tags, 
#    sorted by tag count.

my @Ntagged_all_refs 
    = map { [ $_, $genes_tagged{$_} ] } @Ntagged_genes; 

# N.B.: 'push' alters its first argument!
push ( @Ntagged_all_refs, 
       ( map { [ $_, $ambivs_tagged{$_} ] } @Ntagged_ambivs ), );

push ( @Ntagged_all_refs, 
       ( map { [ $_, $genes_spliced{$_} ] } @unitagged_splices ), );

# See "Note on sorting" above.

my @Ntagged_all_sorted_refs 
    = sort { int($b->[1]) <=> int($a->[1]); } 
      sort { $a->[0] cmp $b->[0]; }
      @Ntagged_all_refs;

foreach my $generef ( @Ntagged_all_sorted_refs ) { 
    print $generef->[0], "\t", $generef->[1], "\n";
}

