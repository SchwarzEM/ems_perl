#!/usr/bin/env perl

use strict;
use warnings;
use List::MoreUtils qw(uniq);

my $family    = $ARGV[0];
my $expr_data = $ARGV[1];

my $data_ref;

open my $EXPR, '<', $expr_data or die "Can't open expression data $expr_data: $!";
while (my $input = <$EXPR>) { 
    chomp $input;
    if ( $input =~ /\A (WBGene\d+) \t (?: [^\t]* \t){5} ([^\t]+) \t (\S+) /xms ) { 
        my $wbgene     = $1;
        my $anatomy    = $2;
        my $anatomy_id = $3;
        my $anatomy_term = "$anatomy [$anatomy_id]";
        $data_ref->{'wb_gene'}->{$wbgene}->{'anatomy_term'}->{$anatomy_term} = 1;
    }
}
close $EXPR or die "Can't close filehandle to expression data $expr_data: $!";

open my $FAM, '<', $family or die "Can't open family table $family: $!";
while (my $input = <$FAM>) {
    chomp $input;
    my $wb_gene_anat_text = q{};
    my @wb_gene_anats = ();
    while ( $input =~ /(WBGene\d+)/xmsg ) { 
        my $wbgene = $1;
        my @more_wb_gene_anats = ();
        @more_wb_gene_anats = sort keys %{ $data_ref->{'wb_gene'}->{$wbgene}->{'anatomy_term'} };
        push @wb_gene_anats, @more_wb_gene_anats;
    }
    @wb_gene_anats = sort @wb_gene_anats;
    @wb_gene_anats = uniq @wb_gene_anats;
    if (@wb_gene_anats) { 
        $wb_gene_anat_text = join '; ', @wb_gene_anats;
    }
    print "$input\t$wb_gene_anat_text\n";
}
close $FAM or die "Can't close filehandle to family table $family: $!";
