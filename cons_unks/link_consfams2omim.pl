#!/usr/bin/env perl

use strict;
use warnings;
use List::MoreUtils qw(uniq);

my $family = $ARGV[0];
my $omim   = $ARGV[1];

my $data_ref;

open my $OMIM, '<', $omim or die "Can't open OMIM genemap file $omim: $!";
while (my $input = <$OMIM>) { 
    chomp $input;

# For the file genemap, the fields are, in order :
# 6  - Gene Symbol(s)
# 10 - MIM Number
# 14 - Disorders (each disorder is followed by its MIM number, if
#         different from that of the locus, and phenotype mapping method (see
#         below).  Allelic disorders are separated by a semi-colon.
# E.g.:
# 1.7|10|2|07|1p35.2|SDC3, SYND3, SDCN|C|Syndecan 3||186357|REc, R, H|||{Obesity, association with}, 601665 (3)| | |4(Synd3)|
# 1.3|12|22|87|1pter-p36|ERPL1, HLM2|C|Endogenous retroviral pol gene-like sequence 1 (oncogene HLM2)||131190|REa, F|||| | ||

    if ( $input =~ /\A (?: [^\|]* \|){5} ([^,\s\|]+) [^\|]*? \| (?: [^\|]* \|){3} (\d+) \| (?: [^\|]* \|){3} ([^\|]*) \| /xms ) { 
        my $hgnc_gene   = $1;
        my $omim_number = $2;
        my $disorders   = $3;
        if ( $disorders =~ / \( 3 \) /xms ) { 
            my $disease = 'OMIM:' . $omim_number . q{|} . $hgnc_gene . q{ [} . $disorders . q{]};
            $data_ref->{'hgnc_gene'}->{$hgnc_gene}->{'disease'}->{$disease} = 1;
        }
    }
    else { 
        die "From OMIM genemap file $omim, can't parse input: $input\n";
    }
}
close $OMIM or die "Can't close filehandle to OMIM genemap file $omim: $!";

open my $FAM, '<', $family or die "Can't open family table $family: $!";
while (my $input = <$FAM>) {
    chomp $input;
    my $hgnc_gene_omim_text = q{};
    my @hgnc_gene_diseases  = ();

    # Sample input:
    # PF15189|DUF4582 13.0013 3 genes c_elegans|WBGene00012648|Y39A1A.9; human|ENSG00000180336|C17orf104; mouse|MGI:2686410|Gm1564 

    while ( $input =~ / human \| ENSG\d+ \| ([^;]+) ; /xmsg ) {
        my $hgnc_gene = $1;
        my @more_hgnc_gene_diseases = ();
        @more_hgnc_gene_diseases = sort keys %{ $data_ref->{'hgnc_gene'}->{$hgnc_gene}->{'disease'} }; 
        push @hgnc_gene_diseases, @more_hgnc_gene_diseases;
    }

    @hgnc_gene_diseases = sort @hgnc_gene_diseases ;
    @hgnc_gene_diseases = uniq @hgnc_gene_diseases ;

    if (@hgnc_gene_diseases) {
        $hgnc_gene_omim_text = join '; ', @hgnc_gene_diseases;
    }

    print "$input\t$hgnc_gene_omim_text\n";
}
close $FAM or die "Can't close filehandle to family table $family: $!";

