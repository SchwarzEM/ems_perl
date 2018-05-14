#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

my $canonical = $ARGV[0];
my $big_table = $ARGV[1];

my $data_ref;

# We will pass this over from the big table to our new output.
my $header = q{};

open my $CANON, '<', $canonical;
while (my $input = <$CANON>) { 
    chomp $input;
    if ( $input =~ /\A (WBGene\d+) \t .+ \t (WBVar\d+) \t .+ \z/xms ) {
        my $gene   = $1;
        my $allele = $2;
        $data_ref->{'canonical_gene'}->{$gene}->{'canonical_allele'} = $allele;
    }
    elsif ( $input !~ /\A Gene \t /xms ) {
        die "Cannot parse input: $input\n";
    }
}
close $CANON;

open my $BIG, '<', $big_table;
while (my $input = <$BIG>) {
    chomp $input;
    if ( $input =~ /\A (WBGene\d+) \t (WBVar\d+) \t .+ \t (WBPhenotype:\d+) \t .+ \z/xms ) { 
        my $gene   = $1;
        my $allele = $2;
        my $pheno  = $3;
        # For each combination of gene, allele, and phenotype, record the full data line.
        $data_ref->{'gene'}->{$gene}->{'allele'}->{$allele}->{'pheno'}->{$pheno} = $input;
    }
    elsif ( $input !~ /\A Gene \t /xms ) {
        die "Cannot parse input: $input\n";
    }
    elsif ( $input =~ /\A Gene \t /xms ) { 
        $header = $input;
    }
}
close $BIG;

my @genes = sort keys %{ $data_ref->{'gene'} };
foreach my $gene (@genes) {
    my $choice_allele = q{};

    if ( exists $data_ref->{'canonical_gene'}->{$gene}->{'canonical_allele'} ) {
        my $possible_choice_allele = $data_ref->{'canonical_gene'}->{$gene}->{'canonical_allele'};
        if ( exists $data_ref->{'gene'}->{$gene}->{'allele'}->{$possible_choice_allele}->{'pheno'} ) { 
            $choice_allele = $possible_choice_allele;
        }
    }

    # if we haven't already picked a $choice allele, then...
    if (! $choice_allele) {
        my %allele2pheno_count = ();
        my @alleles = sort keys %{ $data_ref->{'gene'}->{$gene}->{'allele'} };
        foreach my $allele (@alleles) {
            my @phenos = sort keys %{ $data_ref->{'gene'}->{$gene}->{'allele'}->{$allele}->{'pheno'} };
            my $pheno_count = @phenos;
            $allele2pheno_count{$allele} = $pheno_count;
        }
        my @sorted_alleles = ();
        @sorted_alleles    = sort { $allele2pheno_count{$b} <=> $allele2pheno_count{$a} } @alleles;
        $choice_allele     = $sorted_alleles[0];
    }
    
    if (! $choice_allele) {
        die "Failed to define allele of choice for gene $gene\n";
    }

    my @final_phenos = sort keys %{ $data_ref->{'gene'}->{$gene}->{'allele'}->{$choice_allele}->{'pheno'} };
    my @data_lines = ();
    foreach my $final_pheno (@final_phenos) {
        my $data_line = $data_ref->{'gene'}->{$gene}->{'allele'}->{$choice_allele}->{'pheno'}->{$final_pheno};
        push @data_lines, $data_line;
    }
    @data_lines = sort @data_lines;
    @data_lines = uniq(@data_lines);
    foreach my $data_line (@data_lines) {
        # Print the header only once, at the start of the overall output.
        print "$header\n" if $header;
        $header = q{};

        print "$data_line\n";
    }
}

