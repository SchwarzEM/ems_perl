#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use List::MoreUtils qw(uniq);

my $data_ref;

my $washu_tx2gene = q{} ;
my $washu2acey  = q{};

# Examples of input data files:
# annots/E20201108-05.genes.txt and annots/E20201108-07.genes.txt

my $data_file1 = q{};
my $data_file2 = q{};

$washu_tx2gene = $ARGV[0] if $ARGV[0];
$washu2acey    = $ARGV[1] if $ARGV[1];
$data_file1    = $ARGV[2] if $ARGV[2];
$data_file2    = $ARGV[3] if $ARGV[3];

my $header = "Gene\tES_component\tES_observations\tWashU-spec_ES\tGenes_via_WashU";

my $data_tag1 = $data_file1;
my $data_tag2 = $data_file2;

$data_tag1 =~ s/\Aannots\///;
$data_tag1 =~ s/\.genes\.txt\z//;

$data_tag2 =~ s/\Aannots\///;
$data_tag2 =~ s/\.genes\.txt\z//;

open my $WASHU_TX2GENE, '<', $washu_tx2gene;
while (my $tx2gene = <$WASHU_TX2GENE> ) {
    chomp $tx2gene;
    if ( $tx2gene =~ /\A (\S+) \t (\S+) \z/xms ) {
        my $washu_tx   = $1;
        my $washu_gene = $2;
        $data_ref->{'washu_tx'}->{$washu_tx}->{'washu_gene'} = $washu_gene;    
    }
    else {
        die "From washu_tx2gene input file $washu_tx2gene, cannot parse: $tx2gene\n";
    }
}
close $WASHU_TX2GENE;

open my $WASHU2ACEY, '<', $washu2acey;
while (my $mapping = <$WASHU2ACEY>) {
    chomp $mapping;
    if ( $mapping =~ /\A (\S+) \t (\S+) \z/xms ) {
        my $washu_tx   = $1;
        my $acey_gene  = $2;
        my $washu_gene = q{};

        if ( exists $data_ref->{'washu_tx'}->{$washu_tx}->{'washu_gene'} ) {
            $washu_gene = $data_ref->{'washu_tx'}->{$washu_tx}->{'washu_gene'};
        }
        else {
            die "Cannot map WashU tx $washu_tx to WashU gene\n";
        }

        $data_ref->{'washu_gene'}->{$washu_gene}->{'acey_gene'}->{$acey_gene} = 1;
        $data_ref->{'acey_gene'}->{$acey_gene}->{'washu_gene'}->{$washu_gene} = 1;
    }
    else {
        die "From WashU-to-Acey mapping file $washu2acey, cannot parse: $mapping\n";
    }
}
close $WASHU2ACEY;

open my $DATA_FILE1, '<', $data_file1;
while (my $orig_gene = <$DATA_FILE1>) {
    chomp $orig_gene;

    my @genes = ();
    if ( exists $data_ref->{'washu_gene'}->{$orig_gene} ) {
        @genes = sort keys %{ $data_ref->{'washu_gene'}->{$orig_gene}->{'acey_gene'} };
    }
    else {
        push @genes, $orig_gene;
    }

    foreach my $gene (@genes) {
        $data_ref->{'gene'}->{$gene}->{'data_tag'}->{$data_tag1} = 1;
    }
}
close $DATA_FILE1;

open my $DATA_FILE2, '<', $data_file2;
while (my $orig_gene = <$DATA_FILE2>) {
    chomp $orig_gene;
    my @genes = ();

    if ( exists $data_ref->{'washu_gene'}->{$orig_gene}->{'acey_gene'} ) {
        @genes = sort keys %{ $data_ref->{'washu_gene'}->{$orig_gene}->{'acey_gene'} };
    }
    else {
        push @genes, $orig_gene;
    }

    foreach my $gene (@genes) {
        $data_ref->{'gene'}->{$gene}->{'data_tag'}->{$data_tag2} = 1;
    }
}
close $DATA_FILE2;

my @genes = sort keys %{ $data_ref->{'gene'} };
foreach my $gene (@genes) {
    print "$header\n" if $header;
    $header = q{};

    my @washu_genes = ();
    if ( exists $data_ref->{'acey_gene'}->{$gene}->{'washu_gene'} ) {
        @washu_genes = sort keys %{ $data_ref->{'acey_gene'}->{$gene}->{'washu_gene'} };
    }
    my $WashU_spec_ES = join '; ', @washu_genes;

    my @washu2acey_genes = ();
    foreach my $washu_gene (@washu_genes) {
        if ( exists $data_ref->{'washu_gene'}->{$washu_gene}->{'acey_gene'} ) {
            my @tmp = ();
            @tmp    = sort keys %{ $data_ref->{'washu_gene'}->{$washu_gene}->{'acey_gene'} };
            push @washu2acey_genes, @tmp;
        }
    }
    @washu2acey_genes = sort(@washu2acey_genes);
    @washu2acey_genes =	uniq(@washu2acey_genes);
    my $washu2acey_gene_text = join '; ', @washu2acey_genes;

    my @data_tags = sort keys %{ $data_ref->{'gene'}->{$gene}->{'data_tag'} };
    my $data_tag_text = join '; ', @data_tags;

    print "$gene\tES_component\t$data_tag_text\t$WashU_spec_ES\t$washu2acey_gene_text\n";
}

