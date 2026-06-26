#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $pangene_list   = q{};
my $pg_member_list = q{};
my $taxon2fasta    = q{};
my $suffix         = 'fa';

$pangene_list   = $ARGV[0] if $ARGV[0];
$pg_member_list = $ARGV[1] if $ARGV[1];
$taxon2fasta    = $ARGV[2] if $ARGV[2];
$suffix         = $ARGV[3] if $ARGV[3];

my $data_ref;

if ( (! $pangene_list ) or (! $pg_member_list ) or (! $taxon2fasta ) ) {
    die "Format: extract_syntelogs_26jun2026a.pl [pangene list] [pangene member list] [taxon-to-FASTA table] [optional out-FASTA suffix; default 'fa'] > [extraction line commands]\n";
}

open my $TAXON2FASTA, '<', $taxon2fasta;
while ( my $input = <$TAXON2FASTA> ) {
    chomp $input;
    # Sample input:
    # Aroian    /ocean/projects/mcb190015p/schwarze/necator/2026.01.14/genespace_05/peptide/Aroian.fa
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
        my $taxon = $1;
        my $fasta = $2;
        $data_ref->{'taxon'}->{$taxon}->{'fasta'} = $fasta;
    }
    else {
        die "In taxon2fasta file $taxon2fasta, cannot parse: $input\n";
    }
}
close $TAXON2FASTA;

open my $PG_LIST, '<', $pangene_list;
while ( my $input = <$PG_LIST> ) {
    chomp $input;
    # Sample input:
    # 7961    Aroian; Baylor; Ilik2; Keiser; Mag3; Oita; obscurus     pangene_10030; pangene_10034; pangene_10036;
    if ( $input =~ /\A \S+ \t [^\t]+ \t ([^\t]+) \z/xms ) {
        my $pangene_text = $1;
        my @pangenes = split '; ', $pangene_text;
        foreach my $pangene (@pangenes) {
            $data_ref->{'listed_pg'}->{$pangene} = 1;
        }
    }
    else {
        die "From pangene list $pangene_list, cannot parse: $input\n";
    }
}
close $PG_LIST;

open my $PG_MEMBERS, '<', $pg_member_list;
while ( my $input = <$PG_MEMBERS> ) {
    chomp $input;
    # Sample input:
    # pangene_10030   Aroian|Necator_chrII.g4608; Baylor|Anhui_chrII.g6497; Ilik2|Ilik2_chrII.g5003; Keiser|Keiser_chrII.g6759; \
    # Mag3|Mag3_chrII.g3694; Oita|Oita_chrII.g6548; obscurus|N.sp3_chrII.g5274
    # or
    # 
    if ( $input =~ /\A (\S+) \t (.+) \z/xms ) {
        my $pangene      = $1;
        my $members_text = $2;
        my $outfasta = "$pangene.$suffix";
        my $error    = "$pangene.$suffix.err";

        my @members = split '; ', $members_text;
        if ( exists $data_ref->{'listed_pg'}->{$pangene} ) {
            foreach my $member (@members) {
                if ( $member =~ /\A (\S+) \| (\S+) \z/xms ) {
                    my $taxon = $1;
                    my $seqid = $2;
                    my $fasta = q{};

                    if ( exists $data_ref->{'taxon'}->{$taxon}->{'fasta'} ) {
                        $fasta = $data_ref->{'taxon'}->{$taxon}->{'fasta'};
                    }
                    else {
                        die "From pg_member_list $pg_member_list, cannot map taxon $taxon to FASTA in: $input\n";
                    }
                    print "extract_fasta_subset.pl -r $taxon -l $seqid -f $fasta 1>>$outfasta 2>>$error ;\n";
                }
                else {
                    die "From pg_member_list $pg_member_list, cannot parse member $member in: $input\n";
                }
            }
        }
    }
    else {
        die "From pg_member_list $pg_member_list, cannot parse: $input\n";
    }
}
close $PG_MEMBERS;

