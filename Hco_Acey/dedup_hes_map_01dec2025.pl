#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Basename;
use List::MoreUtils qw(uniq);

my $data_ref;

my $hes_prot = q{};
my $hes_map  = q{};

$hes_prot = $ARGV[0] if $ARGV[0];
$hes_map  = $ARGV[1] if $ARGV[1];

my $rev_hes_map = $hes_map;
$rev_hes_map    = basename($rev_hes_map);
$rev_hes_map    =~ s/\.txt\z//;
$rev_hes_map    = "$rev_hes_map.rev.txt";
$rev_hes_map    = safename($rev_hes_map);

my $rev_hes_prot = $hes_prot;
$rev_hes_prot    = basename($rev_hes_prot);
$rev_hes_prot    =~ s/\.fa\z//;
$rev_hes_prot    = "$rev_hes_prot.rev.fa";
$rev_hes_prot    = safename($rev_hes_prot);

my $hes_id  = q{};
my $hes_seq = q{};

# We want to list only biologically important HES IDs:
my @hes_ids = ();

# Collect remappings (many of which will not actually be different from mappings):
my @remaps = ();

if ( (! $hes_prot ) or (! $hes_map ) ) {
     die "Format: dedup_hes_map_01dec2025.pl [HES proteome] [HES mapping] => [revised HES proteome] (and) [revised HES mapping]\n"
}

open my $PROT, '<', $hes_prot;
while ( my $input = <$PROT> ) {
    chomp $input;
    if ( $input =~ /\A [>] (\S+) /xms ) {
        $hes_id = $1;
    }
    elsif ( $input =~ / \S /xms ) {
        $input =~ s/\s//g;
        if (! $hes_id ) {
            die "Cannot map HES seq to HES id\n";
        }
        $data_ref->{'hes_id'}->{$hes_id}->{'hes_seq'} .= $input;
    }
}
close $PROT;

open my $MAP, '<', $hes_map;
while ( my $input = <$MAP> ) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
        my $gene   = $1;
        my $hes_id = $2;
        if (! exists $data_ref->{'hes_id'} ) {
            die "For gene $gene, cannot map HES ID $hes_id to a sequence\n";
        }
        $data_ref->{'hes_id'}->{$hes_id}->{'gene'}->{$gene} = 1;
        push @hes_ids, $hes_id;
    }
}
close $MAP;

foreach my $key_hes_id (@hes_ids) {
    my $key_hes_seq   = $data_ref->{'hes_id'}->{$key_hes_id}->{'hes_seq'} ;
    my @key_hes_genes = sort keys %{ $data_ref->{'hes_id'}->{$key_hes_id}->{'gene'} };
    my $genecount     = @key_hes_genes;
    if ( $genecount == 1 ) {
       foreach my $key_hes_gene (@key_hes_genes) {
           my $remap = "$key_hes_gene\t$key_hes_id";
           push @remaps, $remap;
           $data_ref->{'key_hes_id'}->{$key_hes_id}->{'key_hes_gene'} = $key_hes_gene;
           $data_ref->{'key_hes_id'}->{$key_hes_id}->{'hes_seq'} = $key_hes_seq;
       }
    }
    elsif ( $genecount >= 2 ) {
       my $i = 0;
       foreach my $key_hes_gene (@key_hes_genes) {
           $i++;
           my $split_hes_id = $key_hes_id . '_' . $i;
           my $remap = "$key_hes_gene\t$split_hes_id";
           push @remaps, $remap;
           $data_ref->{'key_hes_id'}->{$split_hes_id}->{'key_hes_gene'} = $key_hes_gene;
           $data_ref->{'key_hes_id'}->{$split_hes_id}->{'hes_seq'} = $key_hes_seq;
       }
    }
    else {
        die "For HES ID $key_hes_id, cannot reliably count HES genes: $genecount\n";
    }
}

# Remove extra copies of recorded mappings:
@remaps = uniq(@remaps);

open my $REVMAP, '>', $rev_hes_map;
foreach my $remap (@remaps) {
    print $REVMAP "$remap\n";
}
close $REVMAP;

open my $REVPROT, '>', $rev_hes_prot;
my @rev_hes_ids = sort keys %{ $data_ref->{'key_hes_id'} };
foreach my $rev_hes_id (@rev_hes_ids) {
    my $rev_hes_gene = $data_ref->{'key_hes_id'}->{$rev_hes_id}->{'key_hes_gene'};
    my $rev_hes_seq  = $data_ref->{'key_hes_id'}->{$rev_hes_id}->{'hes_seq'};
    my @output_lines = unpack("a60" x (length($rev_hes_seq)/60 + 1), $rev_hes_seq);
    print $REVPROT ">$rev_hes_id  gene=$rev_hes_gene\n";
    foreach my $output_line (@output_lines) { 
        if ($output_line =~ /\S/) { 
            print $REVPROT "$output_line\n";
        }
    }
}
close $REVPROT;

sub safename {
    my $_filename = $_[0];
    my $_orig_filename = $_filename;
    if (-e $_orig_filename) {
        my $_suffix1 = 1;
        $_filename = $_filename . ".$_suffix1";
        while (-e $_filename) {
            $_suffix1++;
            $_filename =~ s/\.\d+\z//xms;
            $_filename = $_filename . ".$_suffix1";
        }
    }
    return $_filename;
}

