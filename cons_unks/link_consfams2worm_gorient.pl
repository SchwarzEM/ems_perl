#!/usr/bin/env perl

use strict;
use warnings;
use List::MoreUtils qw(uniq);
use Scalar::Util qw(looks_like_number);

my $family     = $ARGV[0];
my $gen_orient = $ARGV[1];

my $data_ref;

my $wb_gene = q{};
my $gene1   = q{};

my %flag2int_type = ( Yes => 'known',
                      No  => 'predicted', );

open my $GORIENT, '<', $gen_orient or die "Can't open GeneOrienteer prediction file $gen_orient: $!";
while (my $input = <$GORIENT>) { 
    chomp $input;

    # Sample input:
    # User Input: WBGene00007476
    # Gene name: C09E9.1
    # Gene id: WBGene00007476
    # Gene1   Rank    Known-Interaction       Gene2   Score
    # C09E9.1 1       No      sel-2   13.17
    # C09E9.1 2       No      C30C11.4        13.17
    # C09E9.1 3       No      num-1   1.51

    if ( $input =~ /\A User \s+ Input: \s+ (WBGene\d+) /xms ) {
        $wb_gene = $1;
        $gene1 = q{};
    }

    elsif ( $wb_gene and ( $input =~ /\A Gene \s+ name: \s+ (\S+) /xms ) ) {
        $gene1 = $1;
    }

    elsif ( $wb_gene and $gene1 and ( $input =~ /\A $gene1 \s+ \d+ \s+ (Yes|No) \s+ (\S+) \s+ (\S+) \s* \z/xms ) ) { 
        my $known_int = $1;
        my $gene2     = $2;
        my $score     = $3;
        if (! looks_like_number($score) ) { 
            die "From GeneOrienteer prediction file $gen_orient, can't parse number $score in: $input\n";
        }
        if ( $score >= 5.6 ) {
            my $int_type = $flag2int_type{$known_int};
            my $partner = "$gene2 [$int_type|$score]";
            $data_ref->{'wb_gene'}->{$wb_gene}->{'gene_ori_partner'}->{$partner} = 1;
        }
    }
}
close $GORIENT or die "Can't close filehandle to GeneOrienteer prediction file $gen_orient: $!";

open my $FAM, '<', $family or die "Can't open family table $family: $!";
while (my $input = <$FAM>) {
    chomp $input;
    my $wb_gene_gori_text = q{};
    my @wb_gene_goris = ();
    while ( $input =~ /(WBGene\d+)/xmsg ) { 
        my $wb_gene = $1;
        my @more_wb_gene_goris = ();
        @more_wb_gene_goris = sort keys %{ $data_ref->{'wb_gene'}->{$wb_gene}->{'gene_ori_partner'} };
        push @wb_gene_goris, @more_wb_gene_goris;
    }
    @wb_gene_goris = sort @wb_gene_goris;
    @wb_gene_goris = uniq @wb_gene_goris;
    if (@wb_gene_goris) { 
        $wb_gene_gori_text = join '; ', @wb_gene_goris;
        $wb_gene_gori_text = "GOrient: $wb_gene_gori_text";
    }
    print "$input\t$wb_gene_gori_text\n";
}
close $FAM or die "Can't close filehandle to family table $family: $!";

