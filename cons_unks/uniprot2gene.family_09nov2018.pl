#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $genes = q{};
my $fams  = q{};
my $names = q{};

$genes = $ARGV[0] if $ARGV[0];
$fams  = $ARGV[1] if $ARGV[1];
$names = $ARGV[2] if $ARGV[2];

my $data_ref;

if ( (! $genes) or (! $fams) or (! $names) ) {
    die "Format: uniprot2gene.family_09nov2018.pl [uniprot2genes] [uniprot2families] [families2names] > [output]\n";
}

# uniprot2genes:
# Sample input --
# Q27493	elegans|WBGene00008781|rpoa-2

open my $GENES, '<', $genes;
while (my $input = <$GENES>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t ([^\t]+) \z/xms ) { 
        my $uniprot   = $1;
        my $gene_list = $2;
        my @genes     = split /;[ ]+/, $gene_list;
        foreach my $gene (@genes) {
            $data_ref->{'uniprot'}->{$uniprot}->{'gene'}->{$gene} = 1;
            $data_ref->{'gene'}->{$gene}->{'no_family'} = 1;
        }
    }
    else {
        die "In uniprot2genes table $genes, cannot parse: $input\n";
    }
}
close $GENES;

# families2full_names:
# Sample input -- 
# PF00562	pfam|PF00562|RNA_pol_Rpb2_6

open my $NAMES, '<', $names;
while (my $input = <$NAMES>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
        my $family    = $1;
        my $full_name = $2;
       	if ( exists $data_ref->{'family'}->{$family}->{'full_name'} ) {
            die "Redundant full names for $family: \"$full_name\" versus \"$data_ref->{'family'}->{$family}->{'full_name'}\"\n";
        }
        $data_ref->{'family'}->{$family}->{'full_name'} = $full_name;
    }
    else {
        die "In families2names table $names, cannot parse: $input\n";
    }
}
close $NAMES;

# uniprot2families:
# Sample input -- 
# Q27493	PF06883;PF04563;PF04561;PF04565;PF00562;PF04560;

open my	$FAMS, '<', $fams;

while (my $input = <$FAMS>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
        my $uniprot  = $1;
       	my $fam_text = $2;

        if ( exists $data_ref->{'uniprot'}->{$uniprot}->{'gene'} ) {
            $fam_text =~ s/[;]\z//;
            my @families = split /[;]/, $fam_text;
            my @genes = sort keys %{ $data_ref->{'uniprot'}->{$uniprot}->{'gene'} }; 
            foreach my $gene (@genes) {
                delete $data_ref->{'gene'}->{$gene}->{'no_family'};
                foreach my $family (@families) { 
                    my $fam_name = q{};
                    if ( exists $data_ref->{'family'}->{$family}->{'full_name'} ) {
                        $fam_name = $data_ref->{'family'}->{$family}->{'full_name'};
                        $data_ref->{'gene'}->{$gene}->{'families'}->{$fam_name} = 1;
                    }
                    else {
                        warn "Cannot give family $family a full name\n";
                    }
                }
            }
       	}
    }
    else {
        die "In uniprot2families table $fams, cannot parse: $input\n";
    }
}
close $FAMS;

my @genes = sort keys %{ $data_ref->{'gene'} };
foreach my $gene (@genes) {
    my @families = sort keys %{ $data_ref->{'gene'}->{$gene}->{'families'} };
    foreach my $family (@families) {
        print "$gene\t$family\n";
    }
}

