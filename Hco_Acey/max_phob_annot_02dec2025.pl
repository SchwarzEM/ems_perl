#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $data_ref;

my $phobius = q{};
my $maxiso  = q{};

$phobius = $ARGV[0] if $ARGV[0];
$maxiso  = $ARGV[1] if $ARGV[1];

my $header = "Gene\tMax-iso_seq_ID\tMax-iso_Phobius";

if ( (! $phobius ) or (! $maxiso ) ) {
    die "Format: max_phob_annot_02dec2025.pl [short-format Phobius output] [max-iso FASTA with gene=X] > [gene annot: max-iso seqname; its Phobius]\n";
}

open my $PHOBIUS, '<', $phobius;
while (my $input = <$PHOBIUS>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \s .* \s (\S+) \z/xms ) {
        my $seqname    = $1;
        my $phob_annot = $2;
        if ( exists $data_ref->{'seqname'}->{$seqname} ) {
            die "In Phobius annot file $phobius, redundant annotation for seqname $seqname in: $input\n";
        }
        if ( $seqname ne 'SEQENCE' ) {
            $data_ref->{'seqname'}->{$seqname}->{'phob_annot'} = $phob_annot;
        }
    }
    else {
        die "In Phobius annot file $phobius, cannot parse: $input\n";
    }
}
close $PHOBIUS;

open my $MAXISO, '<', $maxiso;
while (my $input = <$MAXISO>) {
    chomp $input;
    if ( $input =~ /\A [>] (\S+) \b .* gene[=] (\S+) /xms ) {
        my $seqname = $1;
        my $gene    = $2;
        if ( exists $data_ref->{'gene'}->{$gene} ) {
            die "In max-iso FASTA file, redundant annotation for gene $gene in: $input\n";
        }
        $data_ref->{'gene'}->{$gene}->{'seqname'} = $seqname;
        my $phob_annot = $data_ref->{'seqname'}->{$seqname}->{'phob_annot'};
        $data_ref->{'gene'}->{$gene}->{'phob_annot'} = $phob_annot;
    }
    elsif ( $input =~ /\A [>] /xms ) {
        die "In max-iso FASTA file, cannot parse header: $input\n";
    }
}
close $MAXISO; 

my @genes = sort keys %{ $data_ref->{'gene'} };
foreach my $gene (@genes) {
    my $maxiso     = $data_ref->{'gene'}->{$gene}->{'seqname'};
    my $phob_annot = $data_ref->{'gene'}->{$gene}->{'phob_annot'};
    if ($header) {
        print "$header\n";
        $header = q{};
    }
    print "$gene\t$maxiso\t$phob_annot\n";
}
