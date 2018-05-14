#!/usr/bin/env perl

# wbook_TFs_to_newfound_TreeFam_19aug2013.pl -- Erich Schwarz <ems394@cornell.edu>, 8/18/2013.
# Given my first TreeFam tables, plus new psi-BLAST extension, parse the latest official WormBook table so that it lists new psi-BLAST-derived TreeFam data -- worm-human-fly or worm-human (worm-fly was not selected for in the psi-BLASTs).

use strict;
use warnings;

my $data_ref;

my $new_worm_treefams      = $ARGV[0];
my $official_treefam_table = $ARGV[1];
my $orig_table_1           = $ARGV[2];

my $treefam_stem = $official_treefam_table;
$treefam_stem    =~ s/txt\z//;

my $worm_human_fly = $treefam_stem . 'NEW.worm_human_fly.txt';
my $worm_human     = $treefam_stem . 'NEW.worm_human.txt';

$worm_human_fly = safename($worm_human_fly);
$worm_human     = safename($worm_human);

open my $WHF, '>', $worm_human_fly or die "Can't open output for worm-human-fly TF file $worm_human_fly: $!";
open my $WHU, '>', $worm_human     or die "Can't open output for worm-human TF file $worm_human: $!";

open my $NEW_TREEFAMS, '<', $new_worm_treefams or die "Can't open new worm-TreeFam table $new_worm_treefams: $!";
while (my $input = <$NEW_TREEFAMS>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (TF\d+) \z/xms ) {
        my $gene       = $1;
        my $treefam_id = $2;
        $data_ref->{'worm_gene'}->{$gene}->{'TreeFamID'}->{$treefam_id} = 1;
    }
    else { 
        die "Can't parse input from worm-TreeFam table $new_worm_treefams: $input\n";
    }
}

close $NEW_TREEFAMS or die "Can't close filehandle to new worm-TreeFam table $new_worm_treefams: $!";

open my $TREEFAM, '<', $official_treefam_table or die "Can't open TreeFam table $official_treefam_table: $!";
while (my $input = <$TREEFAM>) {
    chomp $input;
    if ( $input =~ /\A (TF\d+) \t (\S+) \t (\S+) \z/xms ) { 
        my $treefam_id = $1;
        my $species    = $2;
        my $gene       = $3;
        $data_ref->{'TreeFamID'}->{$treefam_id}->{'species'}->{$species}->{'gene'}->{$gene} = 1;
    }
    else { 
        die "Can't parse input line from TreeFam table $official_treefam_table: $input\n";
    }
}
close $TREEFAM or die "Can't close filehandle to TreeFam table $official_treefam_table: $!";

open my $ORIG_TABLE, '<', $orig_table_1 or die "Can't open original WormBook TF data table $orig_table_1: $!";
my $header = q{};
while (my $input = <$ORIG_TABLE>) {
    chomp $input; 
    if (! $header) { 
        $header = $input . "\tTreeFam";
        print $WHF "$header\n";
        print $WHU "$header\n";
    }
    if ( $input =~ /\A WBGene\d+ \t (\S+) \t /xms ) { 
        my $gene = $1;
        if ( exists $data_ref->{'worm_gene'}->{$gene}->{'TreeFamID'} ) { 
            my @tf_treefams = sort keys %{ $data_ref->{'worm_gene'}->{$gene}->{'TreeFamID'} };
            foreach my $tf_treefam (@tf_treefams) { 
                my $tfam_text = q{};
                my $output    = q{};
                if ( $data_ref->{'TreeFamID'}->{$tf_treefam}->{'species'}->{'homo_sapiens'} ) {
                    $tfam_text = summarize_treefam($tf_treefam);
                    $output = $input . "\t" . $tfam_text;
                    if ( $data_ref->{'TreeFamID'}->{$tf_treefam}->{'species'}->{'drosophila_melanogaster'} ) {
                        print $WHF "$output\n";
                    }
                    else { 
                        print $WHU "$output\n";
                    }
                }
            }
        }
    }
}

close $WHF or die "Can't close filehandle for worm-human-fly TF file $worm_human_fly: $!"; 
close $WHU or die "Can't close filehandle for worm-human TF file $worm_human: $!";

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

