#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils qw(uniq);
use autodie;

my $family  = q{};
my $phobius = q{};

my $data_ref;
my $help;

GetOptions ( 'family|f=s'     => \$family,
             'phobius|p=s{,}' => \$phobius,
             'help|h'         => \$help, 
);

if ( $help or (! $family) or (! $phobius) ) {
    die "Format: link_consfams2wormphobius_15oct2014.pl\n",
        "    --family|-f      [1 table of conserved unknown gene families, one family per line, with WBGene IDs listed]\n",
        "    --phobius|-p     [1 table linking WBGene IDs to Phobius summaries]\n",
        "    --help|-h        [print this message]\n",
        ;
}

open my $PHOBIUS, '<', $phobius;
while (my $input = <$PHOBIUS>) {
    chomp $input;

    # Sample input lines:
    # Gene    Phobius
    # WBGene00000002|F27C8.1|aat-1    TM(12x)
    # WBGene00000012|C50F2.9|abf-1    SigP
    # WBGene00000019|C24F3.5|abt-1    TM(12x); TM(13x)

    if ( $input =~ /\A (WBGene\d+) \S* \t ([^\t]+) /xms ) { 
        my $wbgene        = $1;
        my $phobius_annot = $2;
        if ( $phobius_annot =~ /\S/xms ) {
            my @phobius_annots = split /;\s+/, $phobius_annot;
            push @{ $data_ref->{'wb_gene'}->{$wbgene}->{'phobius_annots'} }, @phobius_annots; 
        }
    }
}
close $PHOBIUS;

open my $FAM, '<', $family;
while (my $input = <$FAM>) {
    chomp $input;
    my $wb_gene_phobius_text  = q{};
    my @wb_gene_phobius_annots = ();

    while ( $input =~ /(WBGene\d+)/xmsg ) { 
        my $wbgene = $1;
        my @more_wb_gene_phobius_annots = ();
        if ( exists $data_ref->{'wb_gene'}->{$wbgene}->{'phobius_annots'} ) { 
            @more_wb_gene_phobius_annots = @{ $data_ref->{'wb_gene'}->{$wbgene}->{'phobius_annots'} };
        }
        push @wb_gene_phobius_annots, @more_wb_gene_phobius_annots;
    }

    @wb_gene_phobius_annots = sort @wb_gene_phobius_annots;
    @wb_gene_phobius_annots = grep { /\S/ } uniq @wb_gene_phobius_annots;

    if (@wb_gene_phobius_annots) {
        $wb_gene_phobius_text = join '; ', @wb_gene_phobius_annots;
        $wb_gene_phobius_text = "Phobius: $wb_gene_phobius_text";
    }
    print "$input\t$wb_gene_phobius_text\n";
}
close $FAM;

