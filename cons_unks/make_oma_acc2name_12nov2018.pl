#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $oma_groups = q{};
my $oma_descs  = q{};

my $data_ref;

$oma_groups = $ARGV[0] if $ARGV[0]; 
$oma_descs  = $ARGV[1] if $ARGV[1];

if ( (! $oma_groups) or (! $oma_descs) ) {
    die "Format: make_oma_acc2name_12nov2018.pl [oma_groups] [oma_descs] > [oma_acc2name]\n";
}

# Sample OMA groups input file:
# # Orthologous groups from OMA release of All.Jun2018
# # This release has 842789 groups covering 9771653 proteins from 2167 species
# # Format: group number<tab>Fingerprint<tab>tab-separated list of OMA Entry IDs
# 1	AVIWVDE DESA100634	DESM000491	STAHD00638	STAMF00150	THEAM00507	THEGJ01391
# 2	YIRQENI HALLT03360	NATGS03323

open my $OMA_GROUPS, '<', $oma_groups;
while (my $input = <$OMA_GROUPS>) {
    chomp $input;
    if  ( ( $input !~ /\A [#]/xms ) and ( $input =~ /\A (\d+) \t (\S+) \t \S.+\S \z/xms ) ) {
        my $oma_no  = $1;
        my $oma_sig = $2;
        # Some OMA groups have 'n/a' as a signature; these are not usable with UniProt, which maps to defined signatures only.
        if ( $oma_sig ne q{n/a} ) {
            if ( exists $data_ref->{'oma_sig'}->{$oma_sig}->{'oma_no'} ) {
                die "Redundant mapping of OMA signature $oma_sig to OMA numbers $data_ref->{'oma_sig'}->{$oma_sig}->{'oma_no'} and $oma_no\n";
            }
            $data_ref->{'oma_sig'}->{$oma_sig}->{'oma_no'} = $oma_no;
        }
    }
    elsif ( $input =~ /\A [#]/xms ) {
        die "From OMA groups file $oma_groups, cannot parse: $input\n";
    }
}
close $OMA_GROUPS;

# Sample OMA descriptions input file:
# 1	Membrane-anchored protein predicted to be involved in regulation of amylopullulanase-like protein
# 2	N-6 DNA methylase
# 3	Putative uncharacterized protein

open my $OMA_DESCS, '<', $oma_descs;
while (my $input = <$OMA_DESCS>) {
    chomp $input;
    if ( $input =~ /\A (\d+) \t ([^\t]+) \z/xms ) {
        my $oma_no   = $1;
        my $oma_desc = $2;
        $oma_desc =~ s/ /_/g;
        if ( $oma_desc eq q{-} ) {
            $oma_desc = 'unnamed_in_OMA';
        }
        if ( exists $data_ref->{'oma_no'}->{$oma_no}->{'oma_desc'} ) {
            die "Redundant mapping of OMA number $oma_no to two descriptions: \"$data_ref->{'oma_no'}->{$oma_no}->{'oma_desc'}\" and \"$oma_desc\"\n";
        }
        $data_ref->{'oma_no'}->{$oma_no}->{'oma_desc'} = $oma_desc;
    }
    else {
        die "From OMA descriptions file $oma_descs, cannot parse: $input\n";
    }
}
close $OMA_DESCS;

my @oma_sigs = sort keys %{ $data_ref->{'oma_sig'} };
foreach my $oma_sig (@oma_sigs) {
    if (! exists $data_ref->{'oma_sig'}->{$oma_sig}->{'oma_no'} ) {
        die "Cannot map OMA signature $oma_sig to OMA number\n";
    }
    my $oma_no = $data_ref->{'oma_sig'}->{$oma_sig}->{'oma_no'};
    if (! exists $data_ref->{'oma_no'}->{$oma_no}->{'oma_desc'} ) {
        die "Cannot map OMA number $oma_no (and signature $oma_sig) to description\n";
    }
    my $oma_desc = $data_ref->{'oma_no'}->{$oma_no}->{'oma_desc'};
    print "$oma_sig\toma\|$oma_sig\|$oma_desc\n";
}
