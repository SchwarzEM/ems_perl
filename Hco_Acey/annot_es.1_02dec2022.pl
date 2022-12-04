#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

my $data_ref;

# Single-line-per-gene-pair mapping, WashU to v2.1 genes:
my $wa_g2um_g = q{};

# ** E20201108-05 ES genes **

# UMass v2.0 E20201108-05 ES genes:
my $umass_v2_05 = q{};

# UMass v1.0-in-v2.1 E20201108-05 ES genes:
my $umass_v1_05 = q{};

# WashU E20201108-05 ES genes:
my $washu_05 = q{};

# ** E20201108-07 **

# UMass v2.0 E20201108-07 ES genes:
my $umass_v2_07 = q{};

# UMass v1.0-in-v2.1 E20201108-07 ES genes:
my $umass_v1_07 = q{};

# WashU E20201108-07 ES genes:
my $washu_07 = q{};

$wa_g2um_g   = $ARGV[0] if $ARGV[0];

$umass_v2_05 = $ARGV[1] if $ARGV[1];
$umass_v1_05 = $ARGV[2] if $ARGV[2];
$washu_05    = $ARGV[3] if $ARGV[3];

$umass_v2_07 = $ARGV[4] if $ARGV[4];
$umass_v1_07 = $ARGV[5] if $ARGV[5];
$washu_07    = $ARGV[6] if $ARGV[6];

my $header = "Gene\tES_component\tES_observations\tWashU-spec_ES";

if (     ( ! -e $wa_g2um_g ) 
      or ( ! -e $umass_v2_05 ) or ( ! -e $umass_v1_05 ) or ( ! -e $washu_05 ) 
      or ( ! -e $umass_v2_07 ) or ( ! -e $umass_v1_07 ) or ( ! -e $washu_07 ) ) {
    die "Format: annot_es.1_02dec2022.pl ",
        "[WashU to UMass gene table] ",
        "[UMass v2.0 E-05 genes] [UMass v1.0 E-05 genes] [WashU E-05 genes] ",
        "[UMass v2.0 E-07 genes] [UMass v1.0 E-07 genes] [WashU E-07 genes] ",
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

open my $UMASS_V2_05, '<', $umass_v2_05;
while (my $input = <$UMASS_V2_05>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \z/xms ) {
         $data_ref->{'E20201108-05'}->{'umass_gene'}->{$input} = 1;
         $data_ref->{'ES_component'}->{$input} = 1;
    }
    else {
        die "From UMass v2.0 E20201108-05 input file $umass_v2_05, cannot parse: $input\n";
    }
}
close $UMASS_V2_05;

open my $UMASS_V1_05, '<', $umass_v1_05;
while (my $input = <$UMASS_V1_05>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \z/xms ) {
         $data_ref->{'E20201108-05'}->{'umass_gene'}->{$input} = 1;
         $data_ref->{'ES_component'}->{$input} = 1;
    }
    else {
        die "From UMass v1.0 E20201108-05 input file $umass_v1_05, cannot parse: $input\n";
    }
}
close $UMASS_V1_05;

open my $WASHU_05, '<', $washu_05;
while (my $washu_gene = <$WASHU_05>) {
    chomp $washu_gene;
    if ( $washu_gene =~ /\A (\S+) \z/xms ) {
         my @umass_genes = ();
         my $record      = 1;
         if ( exists $data_ref->{'washu_gene'}->{$washu_gene}->{'umass_gene'} ) {
             @umass_genes = sort keys %{ $data_ref->{'washu_gene'}->{$washu_gene}->{'umass_gene'} };
             foreach my $umass_gene (@umass_genes) {
                 if ( exists $data_ref->{'E20201108-05'}->{'umass_gene'}->{$umass_gene} ) {
                     $record = 0;
                 }
             }
             if ($record) {
                 foreach my $umass_gene (@umass_genes) {
                     $data_ref->{'E20201108-05'}->{'umass_via_washu'}->{$umass_gene}->{'washu_gene'}->{$washu_gene} = 1;
                     $data_ref->{'washu_to_umass'}->{$washu_gene}->{'component'}->{'E20201108-05'} = 1;
                     $data_ref->{'ES_component'}->{$umass_gene} = 1;
                 }
             }
         }
         else {
             $data_ref->{'unmappable_washu_gene'}->{$washu_gene}->{'E20201108-05'} = 1
         }
    }
    else {
        die "From WashU E20201108-05 input file $washu_05, cannot parse: $washu_gene\n";
    }
}
close $WASHU_05;

open my $UMASS_V2_07, '<', $umass_v2_07;
while (my $input = <$UMASS_V2_07>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \z/xms ) {
         $data_ref->{'E20201108-07'}->{'umass_gene'}->{$input} = 1;
         $data_ref->{'ES_component'}->{$input} = 1;
    }
    else {
        die "From UMass v2.0 E20201108-07 input file $umass_v2_07, cannot parse: $input\n";
    }
}
close $UMASS_V2_07;

open my $UMASS_V1_07, '<', $umass_v1_07;
while (my $input = <$UMASS_V1_07>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \z/xms ) {
         $data_ref->{'E20201108-07'}->{'umass_gene'}->{$input} = 1;
         $data_ref->{'ES_component'}->{$input} = 1;
    }
    else {
        die "From UMass v1.0 E20201108-07 input file $umass_v1_07, cannot parse: $input\n";
    }
}
close $UMASS_V1_07;

open my $WASHU_07, '<', $washu_07;
while (my $washu_gene = <$WASHU_07>) {
    chomp $washu_gene;
    if ( $washu_gene =~ /\A (\S+) \z/xms ) {
         my @umass_genes = ();
         my $record      = 1;
         if ( exists $data_ref->{'washu_gene'}->{$washu_gene}->{'umass_gene'} ) {
             @umass_genes = sort keys %{ $data_ref->{'washu_gene'}->{$washu_gene}->{'umass_gene'} };
             foreach my $umass_gene (@umass_genes) {
                 if ( exists $data_ref->{'E20201108-07'}->{'umass_gene'}->{$umass_gene} ) {
                     $record = 0;
                 }
             }
             if ($record) {
                 foreach my $umass_gene (@umass_genes) {
                     $data_ref->{'E20201108-07'}->{'umass_via_washu'}->{$umass_gene}->{'washu_gene'}->{$washu_gene} = 1;
                     $data_ref->{'washu_to_umass'}->{$washu_gene}->{'component'}->{'E20201108-07'} = 1;
                     $data_ref->{'ES_component'}->{$umass_gene} = 1;
                 }
             }
         }
         else {
             $data_ref->{'unmappable_washu_gene'}->{$washu_gene}->{'E20201108-07'} = 1;
         }
    }
    else {
        die "From WashU E20201108-07 input file $washu_07, cannot parse: $washu_gene\n";
    }
}
close $WASHU_07;

my @umass_es_genes = sort keys %{ $data_ref->{'ES_component'} };
foreach my $umass_es_gene (@umass_es_genes) {
     my @components = ();
     my @washus     = ();
     if ( exists $data_ref->{'E20201108-05'}->{'umass_gene'}->{$umass_es_gene} ) {
         push @components, 'E20201108-05';
     }
     if ( exists $data_ref->{'E20201108-07'}->{'umass_gene'}->{$umass_es_gene} ) {
         push @components, 'E20201108-07';
     }
     if ( exists $data_ref->{'E20201108-05'}->{'umass_via_washu'}->{$umass_es_gene}->{'washu_gene'} ) {
         push @components, 'E20201108-05';
         my @washus1 = ();
         @washus1    = keys %{ $data_ref->{'E20201108-05'}->{'umass_via_washu'}->{$umass_es_gene}->{'washu_gene'} };
         foreach my $washu1 (@washus1) {
             if (! exists $data_ref->{'washu_to_umass'}->{$washu1}->{'component'} ) {
                 die "Fail to map $washu1 back to its components\n";
             }
             my @components = sort keys %{ $data_ref->{'washu_to_umass'}->{$washu1}->{'component'} };
             my $component_text = join '; ', @components;
             $washu1 = "$washu1 [$component_text]";
             push @washus, $washu1;
         }
     }
     if ( exists $data_ref->{'E20201108-07'}->{'umass_via_washu'}->{$umass_es_gene}->{'washu_gene'} ) {
         push @components, 'E20201108-07';
         my @washus2 = ();
         @washus2    = keys %{ $data_ref->{'E20201108-07'}->{'umass_via_washu'}->{$umass_es_gene}->{'washu_gene'} };
         foreach my $washu2 (@washus2) {
             if (! exists $data_ref->{'washu_to_umass'}->{$washu2}->{'component'} ) {
                 die "Fail to map $washu2 back to its components\n";
             }
             my @components = sort keys %{ $data_ref->{'washu_to_umass'}->{$washu2}->{'component'} };
             my $component_text = join '; ', @components;
             $washu2 = "$washu2 [$component_text]";
             push @washus, $washu2;
         }
     }
     @components = sort(@components);
     @components = uniq(@components);

     @washus = sort(@washus);
     @washus = uniq(@washus);

     my $component_text = q{};
     if (@components) {
         $component_text = join '; ', @components;
     }
     my $washu_text = q{};
     if (@washus) {
         $washu_text = join '; ', @washus;
     }

     print "$header\n" if $header;
     $header = q{};

     print "$umass_es_gene\tES_component\t$component_text\t$washu_text\n";
}

my @unmapped_washu = sort keys %{$data_ref->{'unmappable_washu_gene'} };
foreach my $unmapped_washu (@unmapped_washu) {
    my @components = sort keys %{ $data_ref->{'unmappable_washu_gene'}->{$unmapped_washu} };
    my $component_text = join '; ', @components;
    print "Unmappable\tES_component\t$component_text\t$unmapped_washu\n";
}

