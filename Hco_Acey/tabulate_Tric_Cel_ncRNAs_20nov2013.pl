#!/usr/bin/env perl

use strict;
use warnings;

my @infiles = @ARGV;

my $rna_table = pop @infiles;

my $data_ref;

open my $RNA, '<', $rna_table or die "Can't open RNA table $rna_table: $!";
while (my $input = <$RNA>) { 
    chomp $input;
    if ( $input =~ /\A (\S+) \t (.*) \z/xms ) { 
        my $ncrna = $1;
        my $type  = $2;
        $type     =~ s/\s+\z//;
        $data_ref->{'ncrna'}->{$ncrna}->{'type'}->{$type} = 1;
    }
    else { 
        die "From RNA table $rna_table, can't parse input: $input\n";
    }
}
close $RNA or die "Can't close filehandle to RNA table $rna_table: $!";

foreach my $blast (@infiles) { 
    my $query = q{};
    my $ncrna = q{};
    open my $BLAST, '<', $blast or die "Can't open BlastN report $blast: $!";
    while (my $input = <$BLAST>) {
        chomp $input;
        # Query= scaffold132:626583..626804
        if ( $input =~ /\A Query= [ ] (\S+) /xms ) { 
            $query = $1;
        }
        elsif ( $input =~ /\A > [ ] (\S+) /xms ) {
            $ncrna = $1;
            $data_ref->{'ncrna'}->{$ncrna}->{'tric'}->{$query} = 1;
            $data_ref->{'tric'}->{$query}->{'ncrna'}->{$ncrna} = 1;
        }
    }
    close $BLAST or die "Can't close filehandle to BlastN report $blast: $!";
}

my @tric_dnas = sort keys %{ $data_ref->{'tric'} };
foreach my $tric_dna (@tric_dnas) { 
    my @ncrnas = sort keys %{ $data_ref->{'tric'}->{$tric_dna}->{'ncrna'} };
    foreach my $ncrna (@ncrnas) { 
        my @other_trics = grep { $_ ne $tric_dna } sort keys %{ $data_ref->{'ncrna'}->{$ncrna}->{'tric'} };
        my @ncrna_types = sort keys %{ $data_ref->{'ncrna'}->{$ncrna}->{'type'} };
        my $ncrna_text = join '; ', @ncrna_types;
        $ncrna_text    = "$ncrna [$ncrna_text]";
        my $otric_text = join '; ', @other_trics;
        print "$tric_dna\t$ncrna_text\t$otric_text\n";
    }
}

