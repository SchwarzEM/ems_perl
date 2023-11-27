#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

my $data_ref;

# Single-line-per-gene-pair mapping, WashU v2 to UMass v2.1 genes:
my $wa_g2um_g = q{};

# WashU male ES genes:
my $washu_male = q{};

# WashU female ES genes:
my $washu_female = q{};

$wa_g2um_g    = $ARGV[0] if $ARGV[0];
$washu_male   = $ARGV[1] if $ARGV[1];
$washu_female = $ARGV[2] if $ARGV[2];

my $header = "Gene\tUzoechi_ES";

if (     ( ! -e $wa_g2um_g    ) 
      or ( ! -e $washu_male   )
      or ( ! -e $washu_female ) ) {
    die "Format: annot_es.1_24jan2023.pl ",
        "[WashU to UMass gene table] ",
        "[WashU male ES genes] [WashU female ES genes] ",
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

open my $WASHU_MALE, '<', $washu_male;
while (my $washu_gene = <$WASHU_MALE>) {
    chomp $washu_gene;
    if ( $washu_gene =~ /\A (\S+) \z/xms ) {
         my @umass_genes = ();
         if ( exists $data_ref->{'washu_gene'}->{$washu_gene}->{'umass_gene'} ) {
             @umass_genes = sort keys %{ $data_ref->{'washu_gene'}->{$washu_gene}->{'umass_gene'} };
             foreach my $umass_gene (@umass_genes) {
                 $data_ref->{'male_umass_ES'}->{$umass_gene}->{'washu_gene'}->{$washu_gene} = 1;
                 $data_ref->{'ES_component'}->{$umass_gene} = 1;
             }
         }
         else {
             $data_ref->{'unmappable_washu_gene'}->{$washu_gene}->{'sex'}->{'male'} = 1;
         }
    }
    else {
        die "From WashU male input file $washu_male, cannot parse: $washu_gene\n";
    }
}
close $WASHU_MALE;

open my $WASHU_FEMALE, '<', $washu_female;
while (my $washu_gene = <$WASHU_FEMALE>) {
    chomp $washu_gene;
    if ( $washu_gene =~ /\A (\S+) \z/xms ) {
         my @umass_genes = ();
         if ( exists $data_ref->{'washu_gene'}->{$washu_gene}->{'umass_gene'} ) {
             @umass_genes = sort keys %{ $data_ref->{'washu_gene'}->{$washu_gene}->{'umass_gene'} };
             foreach my $umass_gene (@umass_genes) {
                 $data_ref->{'female_umass_ES'}->{$umass_gene}->{'washu_gene'}->{$washu_gene} = 1;
                 $data_ref->{'ES_component'}->{$umass_gene} = 1;
             }
         }
         else {
             $data_ref->{'unmappable_washu_gene'}->{$washu_gene}->{'sex'}->{'female'} = 1;
         }
    }
    else {
        die "From WashU female ES input file $washu_female, cannot parse: $washu_gene\n";
    }
}
close $WASHU_FEMALE;

my @umass_es_genes = sort keys %{ $data_ref->{'ES_component'} };
foreach my $umass_es_gene (@umass_es_genes) {
     my @washus     = ();

     if ( exists $data_ref->{'male_umass_ES'}->{$umass_es_gene}->{'washu_gene'} ) {
         my @washus1 = ();
         @washus1    = keys %{ $data_ref->{'male_umass_ES'}->{$umass_es_gene}->{'washu_gene'} };
         foreach my $washu1 (@washus1) {
             $washu1 = "$washu1 [Uzoechi_male_ES]";
             push @washus, $washu1;
         }
     }

     if ( exists $data_ref->{'female_umass_ES'}->{$umass_es_gene}->{'washu_gene'} ) {
         my @washus2 = ();
         @washus2    = keys %{ $data_ref->{'female_umass_ES'}->{$umass_es_gene}->{'washu_gene'} };
         foreach my $washu2 (@washus2) {
             $washu2 = "$washu2 [Uzoechi_female_ES]";
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
    my @sexes = sort keys %{ $data_ref->{'unmappable_washu_gene'}->{$unmapped_washu}->{'sex'} };
    my $sex_text = join '; ', @sexes;
    print "Unmappable\tES_component\t$unmapped_washu [$sex_text]\n";
}

