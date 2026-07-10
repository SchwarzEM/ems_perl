#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

my $pg_data  = q{};
my $pg_qvals = q{};

my @pgs = ();

my $header = "Gene\tPangene\tPG_taxa\tp-value\tq-value";

$pg_data  = $ARGV[0] if $ARGV[0];
$pg_qvals = $ARGV[1] if $ARGV[1];

if ( (! $pg_data ) or (! $pg_qvals ) ) {
    die "Format: pangene2annot_09jul2026a.pl [pangene annots] [pangene qvals] > [Aroian genes + taxon combs + qvals]\n";
}

open my $PG_DATA, '<', $pg_data;
while ( my $input = <$PG_DATA> ) {
    chomp $input; 
    # Sample input:
    # pangene_2       Aroian|Necator_chrI.g138; Baylor|Anhui_chrI.g2031; Ilik2|Ilik2_chrI.g365; Keiser|Keiser_chrI.g2488; Oita|Oita_chrI.g2202; obscurus|N.sp3_chrI.g1028
    if ( $input =~ /\A (\S+) \t (\S [^\t]* \S) \z/xms ) {
        my $pg         = $1;
        my $annot_text = $2;
        my $gene       = q{};

        if ( $annot_text =~ / Aroian\|(\S+)/xms ) {
            $gene = $1;
            $gene =~ s/[;]\z//;
        }
        else {
            die "[1] From PG data $pg_data, cannot identify Aroian gene: $input\n";
        }

        if ( exists $data_ref->{'pg'}->{$pg}->{'gene'} ) {
            my $older_name = $data_ref->{'pg'}->{$pg}->{'gene'};
            die "[2] From PG data $pg_data, redundant Aroian genes: $older_name and $gene with $input\n";
        }
        $data_ref->{'pg'}->{$pg}->{'gene'} = $gene;

        while ( $annot_text =~ / (\S+)\| /xmsg ) {
            my $taxon = $1;
            $data_ref->{'pg'}->{$pg}->{'taxon'}->{$taxon} = 1;
        }
    }
    else {
        die "[3] From PG data $pg_data, cannot parse: $input\n";
    }
}
close $PG_DATA;

open my $PG_QVALS, '<', $pg_qvals;
while ( my $input = <$PG_QVALS> ) { 
    chomp $input;

    # Sample inputs:
    # Pangene p-value FDR_BH
    # pangene_12549   0.0     0.0
    # pangene_14877   0.0     0.0

    if ( ( $input =~ /\A (\S+) \t (\S+) \t (\S+) \z/xms ) and ( $input !~ /\A Pangene /xms ) ) {
        my $pg   = $1;
        my $pval = $2;
        my $qval = $3;

        # keep this sorted order for output later
        push @pgs, $pg;

        $data_ref->{'pg'}->{$pg}->{'pval'} = $pval;
        $data_ref->{'pg'}->{$pg}->{'qval'} = $qval;        
    }
    elsif ( $input !~ /\A Pangene /xms ) {
        die "[4] From PG qvals $pg_qvals, cannot parse: $input\n";
    }
}
close $PG_QVALS;

foreach my $pg (@pgs) {
    my $gene = q{};
    my $pval = q{};
    my $qval = q{};
    my @taxa = ();
    my $txtx = q{};    

    if ( exists $data_ref->{'pg'}->{$pg}->{'gene'} ) {
        $gene = $data_ref->{'pg'}->{$pg}->{'gene'};
    }
    if ( exists $data_ref->{'pg'}->{$pg}->{'pval'} ) {
        $pval = $data_ref->{'pg'}->{$pg}->{'pval'};
    }
    if ( exists $data_ref->{'pg'}->{$pg}->{'qval'} ) {
        $qval = $data_ref->{'pg'}->{$pg}->{'qval'};
    }
    if ( exists $data_ref->{'pg'}->{$pg}->{'taxon'} ) {
        @taxa = sort keys %{ $data_ref->{'pg'}->{$pg}->{'taxon'} };
        $txtx = join '; ', @taxa;
    }

    print "$header\n" if $header;
    $header = q{};

    print "$gene\t$pg\t$txtx\t$pval\t$qval\n";
}
