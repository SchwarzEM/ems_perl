#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use Statistics::Descriptive;

my $gene = q{};
my $seq  = q{};

my $data_ref;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A > /xms ) { 
        if ( $input =~ /\A > ((\S+)\.t\d+) /xms ) { 
            $seq  = $1;
            $gene = $2;
            if ( exists $data_ref->{'gene'}->{$gene}->{'seq'}->{$seq} ) { 
                die "Redundant sequence name: $seq\n";
            }
        }
        else { 
            die "Can't parse input line: $input\n";
        }
    }
    elsif ( $input =~ /\S/xms ) { 
        if ( (! $gene ) or (! $seq) ) { 
            die "Cannot map to a named gene and sequence: $input\n";
        }
        $input =~ s/\s//g;
        $data_ref->{'gene'}->{$gene}->{'seq'}->{$seq}->{'residues'} .= $input;
    }
}

my @genes = sort keys %{ $data_ref->{'gene'} };

my $header = "\tMedian_Length\tMedian_GC";
if (@genes) {
    print "$header\n";
}

foreach my $gene1 (@genes) { 
    my @res_counts = ();
    my @gc_fracs   = ();
    my @seqs = sort keys %{ $data_ref->{'gene'}->{$gene1}->{'seq'} };
    foreach my $seq1 (@seqs) {
        if ( exists $data_ref->{'gene'}->{$gene1}->{'seq'}->{$seq1}->{'residues'} ) { 
            my $residues = $data_ref->{'gene'}->{$gene1}->{'seq'}->{$seq1}->{'residues'};
            my $length   = length($residues);
            push @res_counts, $length;

            my $gc_count = ($residues =~ tr/cCgG/cCgG/);
            my $at_count = ($residues =~ tr/aAtT/aAtT/);
            my $nt_count = $gc_count + $at_count;
            my $gc_frac  = ( $gc_count / $nt_count );
            push @gc_fracs, $gc_frac;
        }
        else { 
            die "Failed to identify residues for gene $gene1, sequence $seq1\n";
        }
    }
    my $stat_res_counts = Statistics::Descriptive::Full->new();
    my $stat_gc_fracs = Statistics::Descriptive::Full->new();
    if ( @res_counts and @gc_fracs ) { 
        $stat_res_counts->add_data(@res_counts);
        $stat_gc_fracs->add_data(@gc_fracs);
        my $median_res_count = $stat_res_counts->median();
        my $median_gc_frac   = $stat_gc_fracs->median();
        $median_gc_frac      = sprintf("%.3f", $median_gc_frac);
        print "$gene1\t$median_res_count\t$median_gc_frac\n";
    }
    else { 
        die "Failed to identify residue counts or GC fractions for gene $gene1\n";
    }
}

