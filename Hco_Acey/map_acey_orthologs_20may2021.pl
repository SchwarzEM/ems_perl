#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $cds2gene2species = q{};
my $orthofinder      = q{};
my $immunochip       = q{};

$cds2gene2species = $ARGV[0] if $ARGV[0];
$orthofinder      = $ARGV[1] if $ARGV[1];
$immunochip       = $ARGV[2] if $ARGV[2];

my $front_header = q{};
my $back_header  = q{};

my $data_ref;

if ( (! $cds2gene2species ) or (! $orthofinder ) or (! $immunochip ) ) {
    die "Format: map_acey_orthologs_20may2021.pl [cds2gene2species] [orthofinder] [immunochip] > [orthology-filtered/annotated immunochip] \n";
}

open my $CDS, '<', $cds2gene2species;
while (my $input = <$CDS>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \t (\S+)\z/xms ) {
        my $cds     = $1;
        my $gene    = $2;
        my $taxon   = $3;
        my $genetax = "$gene($taxon)";
        if ( exists $data_ref->{'cds'}->{$cds} ) {
            die "In cds2gene2species file $cds2gene2species, redundant CDS in: $input\n";
        }
        $data_ref->{'cds'}->{$cds}->{'genetax'} = $genetax;
    }
    else {
        die "From cds2gene2species file $cds2gene2species, cannot parse: $input\n";
    }
}
close $CDS;

open my $ICHIP, '<', $immunochip;
while (my $input = <$ICHIP>) {
    chomp $input;

    if ( $input =~ /\A (PROBE_DESIGN_ID \t SEQ_W_POS) \t (.+) \z/xms ) {
        $front_header = $1;
        $back_header  = $2;
    }

    if ( $input =~ /\A (\S+) \t (\S [^\t]* \S) \t .+\z/xms ) {
        my $data_id  = $1;
        my $seq_text = $2;
        $data_ref->{'data_id'}->{$data_id}->{'full_data'} = $input;
        my @seqs = grep { $_ ne 'SEQ_W_POS' } split /;\s+/, $seq_text;
        foreach my $seq (@seqs) {
            my $cds = q{};
            if ( $seq =~ /\A (?: A_duodenale|N_americanus ) _ (\S+):\d+ \z/xms ) {
                $cds = $1;
            }
            else {
                die "In immunochip data file $immunochip, cannot parse hit: $seq\n";
            }

            if (! exists $data_ref->{'cds'}->{$cds} ) {
                die "In immunochip data file $immunochip, for hit $seq, cannot recognize CDS: $cds\n";
            }

            my $genetax = $data_ref->{'cds'}->{$cds}->{'genetax'};
            $data_ref->{'data_id'}->{$data_id}->{'genetax'}->{$genetax} = 1;
            $data_ref->{'genetax'}->{$genetax}->{'data_id'}->{$data_id} = 1;
        }
    }
    else {
	die "In immunochip data file $immunochip, cannot parse: $input\n";
    }
}
close $ICHIP;

open my	$OFIND, '<', $orthofinder;
while (my $input = <$OFIND>) {
    chomp $input;
    if ( $input =~ /\A \S+ \s+ (\S .+ \S)\z/xms ) {
        my $ogroup_txt = $1;
        my @orthologs = split /\s+/, $ogroup_txt;
        my @acey_orths = grep { $_ =~ /\S+\(a_ceylanicum\)/xms } @orthologs;
        my @hkw_orths  = grep {	$_ =~ /\S+\( (?:a_duodenale|necator) \)/xms } @orthologs;
        my $positive   = 0;
        foreach my $hkw_orth (@hkw_orths) {
            if ( exists $data_ref->{'genetax'}->{$hkw_orth}->{'data_id'} ) { 
                foreach my $data_id ( sort keys %{ $data_ref->{'genetax'}->{$hkw_orth}->{'data_id'} } ) {
                    foreach my $acey_orth (@acey_orths) { 
                        $data_ref->{'data_id'}->{$data_id}->{'acey_orth'}->{$acey_orth} = 1;
                        $data_ref->{'data_id_w_acey_orth'}->{$data_id} = 1;
                    }
                }
            }
        }
    }
    else {
        die "In OrthoFinder file $orthofinder, cannot parse: $input\n";
    }
}
close $OFIND;

my @final_data_ids = sort keys %{ $data_ref->{'data_id_w_acey_orth'} };

if (@final_data_ids) {
    print "$front_header\tAcey_homologs\t$back_header\n";
}

foreach my $data_id (@final_data_ids) {
    my $full_data = $data_ref->{'data_id'}->{$data_id}->{'full_data'};
    if ( $full_data =~ /\A (\S+ \t [^\t]+) \t (.+) \z/xms ) { 
        my $front_data = $1;
        my $back_data  = $2;
        my @acey_annots = sort keys %{ $data_ref->{'data_id'}->{$data_id}->{'acey_orth'} }; 
        my $acey_annot_text = join "; ", @acey_annots;
        print "$front_data\t$acey_annot_text\t$back_data\n";
    }
    else {
        die "Cannot parse full data: $full_data\n";
    }
}
