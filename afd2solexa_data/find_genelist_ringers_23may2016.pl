#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Scalar::Util qw(looks_like_number);

my $ref_gene_list = shift @ARGV;

if (! -r $ref_gene_list ) {
    die "Cannot read putative reference gene list: $ref_gene_list\n";
}

my $ref_gene_count = `cat $ref_gene_list | wc -l`;
chomp $ref_gene_count;
if (! looks_like_number($ref_gene_count) ) {
    die "Non-numerical value: \"$ref_gene_count\"\n";
}
$ref_gene_count    = commify($ref_gene_count);

my $header = "Ref. gene list ($ref_gene_list) has $ref_gene_count genes.\n\n"
             . "Indiv_gene_list\tIndiv_gene_count\tOverlap_gene_count\tOverlap_ratio";

foreach my $indiv_gene_list (@ARGV) {
    if (! -r $indiv_gene_list ) { 
        die "Cannot read individual gene list: $indiv_gene_list\n";
    }

    my $indiv_gene_count   = `cat $indiv_gene_list | wc -l`;
    chomp $indiv_gene_count;

    my $overlap_gene_count = `comm -12 $ref_gene_list $indiv_gene_list | wc -l`;
    chomp $overlap_gene_count;

    if (    (! looks_like_number($indiv_gene_count)   ) 
         or (! looks_like_number($overlap_gene_count) ) ) {
        die "Non-numerical value: \"$indiv_gene_count\" or \"$overlap_gene_count\"\n";
    }

    my $overlap_ratio   = 'n/a';
    if ( $indiv_gene_count >= 1 ) {
        $overlap_ratio   = ($overlap_gene_count/$indiv_gene_count);
        $overlap_ratio      = sprintf "%.3f", $overlap_ratio;
    }
    $overlap_gene_count = commify($overlap_gene_count);
    $indiv_gene_count   = commify($indiv_gene_count);
    $overlap_ratio      = "$overlap_ratio [$overlap_gene_count/$indiv_gene_count]";

    print "$header\n" if $header;
    $header = q{};
    print "$indiv_gene_list\t$indiv_gene_count\t$overlap_gene_count\t$overlap_ratio\n";
}

print "\n" if (! $header);


# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}

