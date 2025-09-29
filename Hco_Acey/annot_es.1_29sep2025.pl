#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

my $data_ref;

# Single-line-per-gene-pair mapping, WashU v2 to UMass v2.1 genes:
my $wa_g2um_g = q{};

# WashU L3i ES genes:
my $washu_l3i = q{};

# WashU adult ES genes:
my $washu_adult = q{};

$wa_g2um_g   = $ARGV[0] if $ARGV[0];
$washu_l3i   = $ARGV[1] if $ARGV[1];
$washu_adult = $ARGV[2] if $ARGV[2];

my $header = "Gene\tWong_ES";

if (     ( ! -e $wa_g2um_g    ) 
      or ( ! -e $washu_l3i   )
      or ( ! -e $washu_adult ) ) {
    die "Format: annot_es.1_29sep2025.pl ",
        "[WashU to UMass gene table] ",
        "[WashU L3i ES genes] [WashU adult ES genes] ",
        "> [output]",
        "\n";
}

open my $WA_G2UM_G, '<', $wa_g2um_g;
while (my $input = <$WA_G2UM_G> ) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
        my $washu_gene = $1;
        my $umass_gene = $2;
        $data_ref->{'washu_gene'}->{$washu_gene}->{'umass_gene'}->{$umass_gene} = 1;
    }
    else {
        die "From WashU-to-UMass input file $wa_g2um_g, cannot parse: $input\n";
    }
}
close $WA_G2UM_G;

open my $WASHU_L3I, '<', $washu_l3i;
while (my $washu_gene = <$WASHU_L3I>) {
    chomp $washu_gene;
    if ( $washu_gene =~ /\A (\S+) \z/xms ) {
         my @umass_genes = ();
         if ( exists $data_ref->{'washu_gene'}->{$washu_gene}->{'umass_gene'} ) {
             @umass_genes = sort keys %{ $data_ref->{'washu_gene'}->{$washu_gene}->{'umass_gene'} };
             foreach my $umass_gene (@umass_genes) {
                 $data_ref->{'L3i_umass_ES'}->{$umass_gene}->{'washu_gene'}->{$washu_gene} = 1;
                 $data_ref->{'ES_component'}->{$umass_gene} = 1;
             }
         }
         else {
             $data_ref->{'unmappable_washu_gene'}->{$washu_gene}->{'stage'}->{'L3i'} = 1;
         }
    }
    else {
        die "From WashU male input file $washu_l3i, cannot parse: $washu_gene\n";
    }
}
close $WASHU_L3I;

open my $WASHU_ADULT, '<', $washu_adult;
while (my $washu_gene = <$WASHU_ADULT>) {
    chomp $washu_gene;
    if ( $washu_gene =~ /\A (\S+) \z/xms ) {
         my @umass_genes = ();
         if ( exists $data_ref->{'washu_gene'}->{$washu_gene}->{'umass_gene'} ) {
             @umass_genes = sort keys %{ $data_ref->{'washu_gene'}->{$washu_gene}->{'umass_gene'} };
             foreach my $umass_gene (@umass_genes) {
                 $data_ref->{'adult_umass_ES'}->{$umass_gene}->{'washu_gene'}->{$washu_gene} = 1;
                 $data_ref->{'ES_component'}->{$umass_gene} = 1;
             }
         }
         else {
             $data_ref->{'unmappable_washu_gene'}->{$washu_gene}->{'stage'}->{'adult'} = 1;
         }
    }
    else {
        die "From WashU female ES input file $washu_adult, cannot parse: $washu_gene\n";
    }
}
close $WASHU_ADULT;

my @umass_es_genes = sort keys %{ $data_ref->{'ES_component'} };
foreach my $umass_es_gene (@umass_es_genes) {
     my @washus     = ();

     if ( exists $data_ref->{'L3i_umass_ES'}->{$umass_es_gene}->{'washu_gene'} ) {
         my @washus1 = ();
         @washus1    = keys %{ $data_ref->{'L3i_umass_ES'}->{$umass_es_gene}->{'washu_gene'} };
         foreach my $washu1 (@washus1) {
             $washu1 = "$washu1 [Wong_L3i_ES]";
             push @washus, $washu1;
         }
     }

     if ( exists $data_ref->{'adult_umass_ES'}->{$umass_es_gene}->{'washu_gene'} ) {
         my @washus2 = ();
         @washus2    = keys %{ $data_ref->{'adult_umass_ES'}->{$umass_es_gene}->{'washu_gene'} };
         foreach my $washu2 (@washus2) {
             $washu2 = "$washu2 [Wong_adult_ES]";
             push @washus, $washu2;
         }
     }

     @washus = sort(@washus);
     @washus = uniq(@washus);

     my $washu_text = q{};
     if (@washus) {
         $washu_text = join '; ', @washus;
     }

     print "$header\n" if $header;
     $header = q{};

     print "$umass_es_gene\t$washu_text\n";
}

my @unmapped_washu = sort keys %{ $data_ref->{'unmappable_washu_gene'} };
foreach my $unmapped_washu (@unmapped_washu) {
    my @stages = sort keys %{ $data_ref->{'unmappable_washu_gene'}->{$unmapped_washu}->{'stage'} };
    my $stage_text = join '; ', @stages;
    print "Unmappable\tES_component\t$unmapped_washu [$stage_text]\n";
}

