#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils qw(uniq);
use autodie;

my $family    = q{};
my @express   = ();
my @phenos    = ();
my @no_phenos = ();

my $data_ref;

my $help;

GetOptions ( 'family|f=s'        => \$family,
             'express|e=s{,}'    => \@express,
             'phenos|p=s{,}'     => \@phenos,
             'non_phenos|n=s{,}' => \@no_phenos,
             'help|h'            => \$help, 
);

if ( $help or (! $family) or (! @phenos) or (! @no_phenos) ) {
    die "Format: link_consfams2wormpheno_13oct2014.pl\n",
        "    --family|-f      [1 table of conserved unknown gene families, one family per line, with WBGene IDs listed]\n",
        "    --express|-e     [1+ tables linking WBGene IDs to expression patterns]\n",
        "    --phenos|-p      [1+ tables linking WBGene IDs to observed phenotypes (either Variants or RNAis)]\n",
        "    --non_phenos|-n  [1+ tables linking WBGene IDs to observed *non*-phenotypes (either Variants or RNAis)]\n",
        "    --help|-h        [print this message]\n",
        ;
}

foreach my $expr_data (@express) {
    open my $EXPR, '<', $expr_data;
    while (my $input = <$EXPR>) {
        chomp $input;

        # Sample input from cut -f 2-4 gene_annots/ws245_gene2expression_12oct2014.txt | head --lines=3 
        # "WBGene00001386"	"WBbt:0005813"	"body wall musculature"
        # "WBGene00001386"	"WBbt:0005821"	"vulval muscle"
        # "WBGene00001863"	"WBbt:0005813"	"body wall musculature"

        # Get rid of xace's putting everything in "":
        $input =~ s/["]//g;

        if ( $input =~ /\A [^\t]* \t (WBGene\d+) \t ([^\t]+) \t ([^\t]+) \t /xms ) {
            my $wbgene     = $1;
            my $anatomy_id = $2;   
            my $anatomy    = $3;
            my $anatomy_term = "$anatomy [$anatomy_id]";
            $data_ref->{'wb_gene'}->{$wbgene}->{'anatomy_term'}->{$anatomy_term} = 1;
        }
        else { 
            warn "From expression data table $expr_data, can't parse: $input\n";
        }
    }
    close $EXPR;
}

foreach my $pheno (@phenos) {
    open my $PHENO, '<', $pheno;
    while (my $input = <$PHENO>) {
        chomp $input;

        # Sample input lines from cut -f 2,4 gene_annots/ws245_gene2phenos_12oct2014.tsv.txt ;
        # "WBGene00003883"	"amphid phasmid morphology variant"
        # "WBGene00003883"	"male ray morphology variant"
        # "WBGene00003883"	"male mating efficiency reduced"

        # Get rid of xace's putting everything in "":
        $input =~ s/["]//g;

        if ( $input =~ /\A [^\t]* \t (WBGene\d+) \t (?: [^\t]* \t) ([^\t]*) \t /xms ) { 
            my $wbgene    = $1;
            my $phenotype = $2;
            if ( $phenotype =~ /\S/xms ) { 
                $data_ref->{'wb_gene'}->{$wbgene}->{'phenotype'}->{$phenotype} = 1;
            }
        }
        else { 
            warn "From phenotype data $pheno, can't parse: $input\n";
        }
    }
    close $PHENO;
}

foreach my $no_pheno (@no_phenos) {
    open my $NO_PHENO, '<', $no_pheno;
    while (my $input = <$NO_PHENO>) {
        chomp $input;
        
        # Get rid of xace's putting everything in "":
        $input =~ s/["]//g;
            
        if ( $input =~ /\A [^\t]* \t (WBGene\d+) \t (?: [^\t]* \t) ([^\t]*) \t /xms ) {
            my $wbgene    = $1;
            my $no_pheno = $2;
            if ( $no_pheno =~ /\S/xms ) {
                $data_ref->{'wb_gene'}->{$wbgene}->{'no_pheno'}->{$no_pheno} = 1;
            }
        }
        else {
            warn "From non-phenotype data $no_pheno, can't parse: $input\n";
        }
    }
    close $NO_PHENO; 
}

open my $FAM, '<', $family or die "Can't open family table $family: $!";
while (my $input = <$FAM>) {
    chomp $input;

    my $wb_gene_pheno_text    = q{};
    my $wb_gene_no_pheno_text = q{};
    my $wb_gene_anat_text = q{};

    my @wb_gene_phenos = ();
    my @wb_gene_no_phenos = ();
    my @wb_gene_anats = ();

    while ( $input =~ /(WBGene\d+)/xmsg ) { 
        my $wbgene = $1;

        my @more_wb_gene_anats = ();
        my @more_wb_gene_phenos = ();
        my @more_wb_gene_no_phenos = ();

        @more_wb_gene_anats = sort keys %{ $data_ref->{'wb_gene'}->{$wbgene}->{'anatomy_term'} };
        @more_wb_gene_phenos    = sort keys %{ $data_ref->{'wb_gene'}->{$wbgene}->{'phenotype'} };
        @more_wb_gene_no_phenos = sort keys %{ $data_ref->{'wb_gene'}->{$wbgene}->{'no_pheno'} };

        push @wb_gene_anats, @more_wb_gene_anats;
        push @wb_gene_phenos, @more_wb_gene_phenos;
        push @wb_gene_no_phenos, @more_wb_gene_no_phenos;
    }

    @wb_gene_anats = sort @wb_gene_anats;
    @wb_gene_anats = grep { /\S/ } uniq @wb_gene_anats;

    @wb_gene_phenos = sort @wb_gene_phenos;
    @wb_gene_phenos = grep { /\S/ } uniq @wb_gene_phenos;

    @wb_gene_no_phenos = sort @wb_gene_no_phenos;
    @wb_gene_no_phenos = grep { /\S/ } uniq @wb_gene_no_phenos;

    if (@wb_gene_anats) {
        $wb_gene_anat_text = join '; ', @wb_gene_anats;
    }

    if (@wb_gene_phenos) {
        $wb_gene_pheno_text = join '; ', @wb_gene_phenos;
        $wb_gene_pheno_text = "WBPheno: $wb_gene_pheno_text";
    }

    if (@wb_gene_no_phenos) {
        $wb_gene_no_pheno_text = join '; ', @wb_gene_no_phenos;
        $wb_gene_no_pheno_text =  "WBNonPheno: $wb_gene_no_pheno_text";
    }       
    
    print "$input\t$wb_gene_anat_text\t$wb_gene_pheno_text\t$wb_gene_no_pheno_text\n";
}
close $FAM or die "Can't close filehandle to family table $family: $!";

