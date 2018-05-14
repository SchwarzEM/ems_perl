#!/usr/bin/env perl

use strict;
use warnings;

my $x_axis       = $ARGV[0];
my $y_axis       = $ARGV[1];
my $orthologs    = $ARGV[2];
my $housekeeping = $ARGV[3];
my $xenologs     = $ARGV[4];
my $ribosomal    = $ARGV[5];
my $expressed    = $ARGV[6];

my %is_orth   = ();
my %is_house  = ();
my %is_xeno   = ();
my %is_ribo   = ();
my %expressed = ();

my $data_ref;

open my $ORTHS, '<', $orthologs or die "Can't open ortholog list $orthologs: $!";
while (my $input = <$ORTHS>) { 
    chomp $input;
    if ( $input =~ /\A (\S+) \z/xms ) { 
        $is_orth{$input} = 1;
    } 
    else { 
        die "Can't parse input from ortholog list $orthologs: $input\n";
    }
}
close $ORTHS or die "Can't close filehandle to ortholog list $orthologs: $!";

open my $HOUSE, '<', $housekeeping or die "Can't open housekeeping list $housekeeping: $!";
while (my $input = <$HOUSE>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \z/xms ) {
        $is_house{$input} = 1;
    }
    else {
        die "Can't parse input from housekeeping list $housekeeping: $input\n";
    } 
}
close $HOUSE or die "Can't close filehandle to housekeeping list $housekeeping: $!";

open my $XENO, '<', $xenologs or die "Can't open xenologs list $xenologs: $!";
while (my $input = <$XENO>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \z/xms ) {
        $is_xeno{$input} = 1;
    }
    else {
        die "Can't parse input from xenologs list $xenologs: $input\n";
    }
}
close $XENO or die "Can't close filehandle to xenologs list $xenologs: $!";

open my $RIBO, '<', $ribosomal or die "Can't open ribosomal list $ribosomal: $!";
while (my $input = <$RIBO>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \z/xms ) {
        $is_ribo{$input} = 1;
    }
    else {
        die "Can't parse input from ribosomal list $ribosomal: $input\n";
    }
}
close $RIBO or die "Can't close filehandle to ribosomal list $ribosomal: $!";

open my $EXPR, '<', $expressed or die "Can't open expressed list $expressed: $!";
while (my $input = <$EXPR>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \z/xms ) {
        $expressed{$input} = 1;
    }
    else {
        die "Can't parse input from expressed list $expressed: $input\n";
    }
}
close $EXPR or die "Can't close filehandle to xenologs list $xenologs: $!";

open my $X_AXIS, '<', $x_axis or die "Can't open x-axis list $x_axis: $!";
while (my $input = <$X_AXIS>) {
    chomp $input;
    if ( $input =~ /\A a(\S+)\.t\d+ \s+ (\S+) \s* \z/xms ) {
        my $gene    = $1;
        my $x_value = $2;
        $gene = 'Acey_s' . $gene;
        if ( $expressed{$gene} ) { 
            $data_ref->{'gene'}->{$gene}->{'x'} = $x_value;
        }
    }
    else {
        warn "Can't parse input from x-axis list $x_axis: $input\n";
    }
}
close $X_AXIS or die "Can't close filehandle to x-axis list $x_axis: $!";

open my $Y_AXIS, '<', $y_axis or die "Can't open y-axis list $y_axis: $!";
while (my $input = <$Y_AXIS>) {
    chomp $input;
    if ( $input =~ /\A a(\S+)\.t\d+ \s+ (\S+) \s* \z/xms ) {
        my $gene    = $1;
        my $y_value = $2;
        $gene = 'Acey_s' . $gene;
        if ( $expressed{$gene} ) {
            $data_ref->{'gene'}->{$gene}->{'y'} = $y_value;
        }
    }
    else {
        warn "Can't parse input from y-axis list $y_axis: $input\n";
    }
}    
close $Y_AXIS or die "Can't close filehandle to y-axis list $y_axis: $!";

my $header = "\t\tExpressed_gene\tOrtholog\tHousekeeping\tRibosomal\tXenolog";

my @genes = sort keys %{ $data_ref->{'gene'} };

foreach my $gene (@genes) { 
    my $x_value = $data_ref->{'gene'}->{$gene}->{'x'};
    my $y_value = $data_ref->{'gene'}->{$gene}->{'y'};
    my $i = 1;
    if ($is_orth{$gene}) {
        $i = 2;
    }
    if ($is_house{$gene}) {
        $i = 3;
    }
    if ($is_ribo{$gene}) {
        $i = 4;
    }
    if ($is_xeno{$gene}) {
        $i = 5;
    }
    if ($header) { 
        print "$header\n";
        $header = q{};
    }
    my $spacer = ( "\t" x $i );
    print "$gene\t$x_value$spacer$y_value\n";
}

