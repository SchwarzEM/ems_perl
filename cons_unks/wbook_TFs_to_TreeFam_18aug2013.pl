#!/usr/bin/env perl

# wbook_TFs_to_TreeFam_18aug2013.pl -- Erich Schwarz <ems394@cornell.edu>, 8/18/2013.
# Given my first TreeFam tables, parse the latest official WormBook table so that it has TreeFam data where useful -- worm-fly or worm-human, etc.

use strict;
use warnings;

my $data_ref;

my $treefam_table = $ARGV[0];
my $orig_table_1  = $ARGV[1];

my $treefam_stem = $treefam_table;
$treefam_stem    =~ s/txt\z//;

my $worm_human_fly = $treefam_stem . 'worm_human_fly.txt';
my $worm_human     = $treefam_stem . 'worm_human.txt';
my $worm_fly       = $treefam_stem . 'worm_fly.txt';

$worm_human_fly = safename($worm_human_fly);
$worm_human     = safename($worm_human);
$worm_fly       = safename($worm_fly);

open my $WHF, '>', $worm_human_fly or die "Can't open output for worm-human-fly TF file $worm_human_fly: $!";
open my $WHU, '>', $worm_human     or die "Can't open output for worm-human TF file $worm_human: $!";
open my $WFY, '>', $worm_fly       or die "Can't open output for worm-fly TF file $worm_fly: $!";

open my $TREEFAM, '<', $treefam_table or die "Can't open TreeFam table $treefam_table: $!";
while (my $input = <$TREEFAM>) {
    chomp $input;
    if ( $input =~ /\A (TF\d+) \t (\S+) \t (\S+) \z/xms ) { 
        my $treefam_id = $1;
        my $species    = $2;
        my $gene       = $3;
        $data_ref->{'TreeFamID'}->{$treefam_id}->{'species'}->{$species}->{'gene'}->{$gene} = 1;
        if ( $species eq 'caenorhabditis_elegans' ) { 
            $data_ref->{'worm_gene'}->{$gene}->{'TreeFamID'}->{$treefam_id} = 1;
        }
    }
    else { 
        die "Can't parse input line from TreeFam table $treefam_table: $input\n";
    }
}
close $TREEFAM or die "Can't close filehandle to TreeFam table $treefam_table: $!";

open my $ORIG_TABLE, '<', $orig_table_1 or die "Can't open original WormBook TF data table $orig_table_1: $!";
my $header = q{};
while (my $input = <$ORIG_TABLE>) {
    chomp $input; 
    if (! $header) { 
        $header = $input . "\tTreeFam";
        print $WHF "$header\n";
        print $WHU "$header\n";
        print $WFY "$header\n";
    }
    if ( $input =~ /\A WBGene\d+ \t (\S+) \t /xms ) { 
        my $gene = $1;
        if ( exists $data_ref->{'worm_gene'}->{$gene}->{'TreeFamID'} ) { 
            my @tf_treefams = sort keys %{ $data_ref->{'worm_gene'}->{$gene}->{'TreeFamID'} };
            foreach my $tf_treefam (@tf_treefams) { 
                my $tfam_text = q{};
                my $output    = q{};
                if (   ( $data_ref->{'TreeFamID'}->{$tf_treefam}->{'species'}->{'homo_sapiens'} )
                    or ( $data_ref->{'TreeFamID'}->{$tf_treefam}->{'species'}->{'drosophila_melanogaster'} ) ) {
                    $tfam_text = summarize_treefam($tf_treefam);
                    $output = $input . "\t" . $tfam_text;
                }
                if (    ( $data_ref->{'TreeFamID'}->{$tf_treefam}->{'species'}->{'homo_sapiens'} ) 
                    and ( $data_ref->{'TreeFamID'}->{$tf_treefam}->{'species'}->{'drosophila_melanogaster'} ) ) {
                    print $WHF "$output\n";
                }
                elsif ( $data_ref->{'TreeFamID'}->{$tf_treefam}->{'species'}->{'homo_sapiens'} ) {
                    print $WHU "$output\n";
                }
                elsif ( $data_ref->{'TreeFamID'}->{$tf_treefam}->{'species'}->{'drosophila_melanogaster'} ) {
                    print $WFY "$output\n";
                }
            }
        }
    }
}

close $WHF or die "Can't close filehandle for worm-human-fly TF file $worm_human_fly: $!"; 
close $WHU or die "Can't close filehandle for worm-human TF file $worm_human: $!";
close $WFY or die "Can't close filehandle for worm-fly TF file $worm_fly: $!";

sub summarize_treefam { 
    my $_treefam_id = $_[0];
    my @_spp_and_genes = ();
    my @_treefam_spp = sort keys %{ $data_ref->{'TreeFamID'}->{$_treefam_id}->{'species'} };
    foreach my $_tfam_species (@_treefam_spp) {
        my @_tf_sp_genes = sort keys %{ $data_ref->{'TreeFamID'}->{$_treefam_id}->{'species'}->{$_tfam_species}->{'gene'} };
        my $_tf_sp_gene_text = join ', ', @_tf_sp_genes;
        $_tf_sp_gene_text = '[' . "$_tfam_species" . '] ' . $_tf_sp_gene_text;
        push @_spp_and_genes, $_tf_sp_gene_text;
    }
    my $_tfam_summary = join '; ', @_spp_and_genes;
    $_tfam_summary = "$_treefam_id" . ": " . $_tfam_summary;
    return $_tfam_summary;
}

sub safename {
    my $filename = $_[0];
    my $orig_filename = $filename;
    if (-e $orig_filename) {
        my $suffix1 = 1;
        $filename = $filename . ".$suffix1";
        while (-e $filename) {
            $suffix1++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix1";
        }
    }
    return $filename;
}

