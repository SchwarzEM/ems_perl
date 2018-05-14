#!/usr/bin/env perl

use strict;
use warnings;

my $asp_full_set = $ARGV[0];
my $asp_aligned  = $ARGV[1];

my $data_ref;

my $header = "Gene\tASP";

my %synonyms = ( 
    'ASP-1'  => 'Acey_2012.08.05_0248.g90',
    'ASP-2'  => 'Acey_2012.08.05_0184.g991',
    'ASP-3A' => 'Acey_2012.08.05_0067.g52',
    'ASP-3B' => 'Acey_2012.08.05_0067.g57',
    'ASP-4A' => 'Acey_2012.08.05_0004.g1750',
    'ASP-4B' => 'Acey_2012.08.05_0004.g1747',
    'ASP-5A' => 'Acey_2012.08.05_0003.g1624',
    'ASP-5B' => 'Acey_2012.08.05_0003.g1625',
    'ASP-6'  => 'Acey_2012.08.05_0067.g36',
    'ASP-7A' => 'Acey_2012.08.05_0222.g2615',
    'ASP-7B' => 'Acey_2012.08.05_0012.g1631',
);

if ( (! $asp_full_set) or (! $asp_aligned) ) { 
    die "Format: build_ASP_set_2013.11.30.pl [full ASP gene list] [aligned ASP gene list]\n";
}

open my $FULL, '<', $asp_full_set or die "Can't open ASP full set, $asp_full_set: $!";
while (my $input = <$FULL>) {
    chomp $input;
    # Sample input;
    # >Acey_2012.08.05_0001.g198.t2
    if ( $input =~ /\A > (\S+) \. t \d+ \b/xms ) { 
        my $asp_gene = $1;
        $data_ref->{'asp'}->{$asp_gene}->{'observed'} = 1;
    }
    elsif ( $input =~ /\A > /xms ) {
        die "From ASP aligned set, $asp_aligned, can't parse: $input\n";
    }
}
close $FULL or die "Can't close filehandle to ASP full set, $asp_full_set: $!";

open my $ALIGN, '<', $asp_aligned or die "Can't open ASP aligned set, $asp_aligned: $!";
while (my $input = <$ALIGN>) { 
    chomp $input;
    # Sample input;
    # >ASP-s0061.g3273/4-187
    if ( $input =~ /\A > (ASP\-\S+) \/ (\d+ \- \d+) \b/xms ) { 
        my $asp_alias = $1;
        my $asp_gene  = q{};
        my $domain    = $2;
        if (! $synonyms{$asp_alias} ) {
            if ( $asp_alias =~ /\A ASP\-s (\d+ \. g \d+) \z/xms ) {
                $asp_gene = $1;
                $asp_gene = 'Acey_2012.08.05_' . $asp_gene;
            }
            else {
                die "From ASP aligned set, $asp_aligned, with alias \"$asp_alias\" and domain \"$domain\", can't parse: $input\n"
            }
        }
        if ( $synonyms{$asp_alias} ) {
            $asp_gene = $synonyms{$asp_alias};
        }
        $data_ref->{'asp'}->{$asp_gene}->{'aligned'}->{$domain} = 1;
        $data_ref->{'asp'}->{$asp_gene}->{'alias'} = $asp_alias;
    }
}
close $ALIGN or die "Can't close filehandle to ASP aligned set, $asp_aligned: $!";

my @asp_genes = sort keys %{ $data_ref->{'asp'} };

foreach my $asp_gene (@asp_genes) {
    # Print header only once, at start:
    if ($header) {
        print "$header\n";
        $header = q{};
    }

    # Default data column text.
    my $data = 'ASP';

    # If this was an aligned ASP:
    if ( exists $data_ref->{'asp'}->{$asp_gene}->{'aligned'} ) {
        my @domains = sort { by_residue($a) <=> by_residue($b) } sort keys %{ $data_ref->{'asp'}->{$asp_gene}->{'aligned'} };
        my $domain_text = join q{, }, @domains;
        my $alias = $data_ref->{'asp'}->{$asp_gene}->{'alias'};
        $data = "ASP [aligned, as \"$alias\", residues $domain_text]";
    }
    print "$asp_gene\t$data\n";
}

sub by_residue {
    my $_block = $_[0];
    if ( $_block =~ /\A (\d+) \- \d+ \z/xms ) { 
        my $_residue = $1;
        return $_residue;
    }
    else { 
        die "Subroutine by_residue failed to sort block $_block\n";
    }
    return;
}
