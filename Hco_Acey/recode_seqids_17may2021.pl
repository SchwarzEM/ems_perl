#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $infile = q{};

$infile = $ARGV[0] if $ARGV[0];

if (! $infile ) {
    die "Format: recode_seqids_17may2021.pl [infile table] > [outfile table with reworked probe data columns]\n";
}

my $data_ref;

open my $INFILE, '<', $infile;
while (my $input = <$INFILE>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \t (\S+) \t (\S+) \t .+ \z/xms ) {
        my $probe_design_id = $1;
        my $seq_id          = $2;
        my $probe_sequence  = $3;
        my $position        = $4;

        if ( ( $seq_id =~ /\AA_duodenale/xms ) or ( $seq_id =~ /\AN_americanus/xms ) ) {
            my $seq_w_pos = "$seq_id:$position";

            if (     ( exists $data_ref->{'probe_design_id'}->{$probe_design_id}->{'probe_sequence'}             ) 
                 and ( $probe_sequence ne $data_ref->{'probe_design_id'}->{$probe_design_id}->{'probe_sequence'} ) 
               ) {
                die "For $probe_design_id, two different probe sequences:",
                    " $data_ref->{'probe_design_id'}->{$probe_design_id}->{'probe_sequence'} and $probe_sequence\n";
            }

            if ( $probe_sequence !~ /\A [ACDEFGHIKLMNPQRSTVWY]+ \z/xms ) {
                die "Aberrant probe sequence: $probe_sequence\n";
            }

            $data_ref->{'probe_design_id'}->{$probe_design_id}->{'probe_sequence'} = $probe_sequence;
            $data_ref->{'probe_design_id'}->{$probe_design_id}->{'seq_w_pos'}->{$seq_w_pos} = 1;            
        }
    }
    else {
        die "When extracting probe data from input file $infile, cannot parse: $input\n";
    }
}
close $INFILE;

open $INFILE, '<', $infile;
while (my $input = <$INFILE>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t \S+ \t \S+ \t \S+ \t (.+) \z/xms ) {
        my $probe_design_id = $1;
        my $bulk_data       = $2;

        if ( $probe_design_id eq 'PROBE_DESIGN_ID' ) {
            print "PROBE_DESIGN_ID\tSEQ_W_POS\tPROBE_SEQUENCE\t$bulk_data\n";
        }
        elsif ( exists $data_ref->{'probe_design_id'}->{$probe_design_id} ) {
            my $probe_sequence = $data_ref->{'probe_design_id'}->{$probe_design_id}->{'probe_sequence'};
            my @seqs_w_pos = sort keys %{ $data_ref->{'probe_design_id'}->{$probe_design_id}->{'seq_w_pos'} };
            my $seq_w_pos_data = join '; ', @seqs_w_pos;
            print "$probe_design_id\t$seq_w_pos_data\t$probe_sequence\t$bulk_data\n";
        }
        else {
            die "Cannot identify probe_design_id \"$probe_design_id\" in input line: $input\n";
        }
    }
    else {
        die "When rereading input file $infile, cannot parse: $input\n";
    }
}
close $INFILE;

