#!/usr/bin/env perl

use strict;
use warnings;

my $aspr_full_set = $ARGV[0];
my $aspr_aligned  = $ARGV[1];
my $aspr_orthomcl = $ARGV[2];

my $data_ref;

my $header = "Gene\tASPR";

if ( (! $aspr_full_set) or (! $aspr_aligned) or (! $aspr_orthomcl) ) { 
    die "Format: build_ASPR_91set_2013.11.30.pl [full ASPR gene list] [aligned ASPR gene list] [ORTHOMCL-member ASPR gene list]\n";
}

open my $FULL, '<', $aspr_full_set or die "Can't open ASPR full set, $aspr_full_set: $!";
while (my $input = <$FULL>) {
    chomp $input;
    # Sample input;
    # >Acey_2012.08.05_0002.g551.t2
    if ( $input =~ /\A > (\S+) \. t \d+ \b/xms ) { 
        my $aspr_gene = $1;
        $data_ref->{'aspr'}->{$aspr_gene}->{'observed'} = 1;
    }
}
close $FULL or die "Can't close filehandle to ASPR full set, $aspr_full_set: $!";

open my $ALIGN, '<', $aspr_aligned or die "Can't open ASPR aligned set, $aspr_aligned: $!";
while (my $input = <$ALIGN>) { 
    chomp $input;
    # Sample input;
    # >ASPR-s0042.g717/7-159
    if ( $input =~ /\A > (ASPR-s (\d+\.g\d+)) \/ (\d+ \- \d+) /xms ) { 
        my $aspr_alias = $1;
        my $aspr_gene  = $2;
        my $domain     = $3;
        $aspr_gene = 'Acey_2012.08.05_' . $aspr_gene;
        $data_ref->{'aspr'}->{$aspr_gene}->{'aligned'}->{$domain} = 1;
        $data_ref->{'aspr'}->{$aspr_gene}->{'alias'} = $aspr_alias;
    }
}
close $ALIGN or die "Can't close filehandle to ASPR aligned set, $aspr_aligned: $!";

open my $ORTHO, '<', $aspr_orthomcl or die "Can't open ASPR ORTHOMCL set, $aspr_orthomcl: $!";
while (my $input = <$ORTHO>) {
    chomp $input;
    # Sample input;
    # Acey_2012.08.05_0015.g2865.t
    if ( $input =~ /\A (\S+) \. t \z/xms ) {
        my $aspr_gene = $1;
        $data_ref->{'aspr'}->{$aspr_gene}->{'orthomcl'} = 1;
    }
}
close $ORTHO or die "Can't close filehandle to ASPR ORTHOMCL set, $aspr_orthomcl: $!";

my @aspr_genes = sort keys %{ $data_ref->{'aspr'} };

foreach my $aspr_gene (@aspr_genes) {
    # Print header only once, at start:
    if ($header) {
        print "$header\n";
        $header = q{};
    }

    # Default data column text.
    my $data = 'ASPR';

    # Two possible options, if this was an aligned ASPR:
    if ( exists $data_ref->{'aspr'}->{$aspr_gene}->{'aligned'} ) {
        my @domains     = sort keys %{ $data_ref->{'aspr'}->{$aspr_gene}->{'aligned'} };
        my $domain_text = join q{, }, @domains;
        my $alias       = $data_ref->{'aspr'}->{$aspr_gene}->{'alias'};
        if ( exists $data_ref->{'aspr'}->{$aspr_gene}->{'orthomcl'} ) {
            $data = "ASPR [original member of ORTHOMCL896.14spp; aligned, as \"$alias\", residues $domain_text]";
        }
        else {
            $data = "ASPR [aligned, as \"$alias\", residues $domain_text]";
        }
    }
    # One remaining option:
    elsif ( exists $data_ref->{'aspr'}->{$aspr_gene}->{'orthomcl'} ) {
        $data = "ASPR [original member of ORTHOMCL896.14spp]";
    }
    print "$aspr_gene\t$data\n";
}

