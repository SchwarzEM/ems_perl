#!/usr/bin/env perl

use strict ;
use warnings;
use autodie;

my $cds2gene = q{};
my $cds_loj  = q{};

$cds2gene = $ARGV[0] if $ARGV[0];
$cds_loj  = $ARGV[1] if $ARGV[1];

my $data_ref;   # declare as completely unspecified variable, *not* as ' = q{}', or data structures will crash under 'use strict'.

# Format of input $cds2gene: "CDS/transcript name \t gene name"
# Extract this from, e.g., ParaSite protein FASTA header lines such as:
# 
# >augustus-ANCCEYDFT_Contig1-processed-gene-17.0-mRNA-1 transcript=augustus-ANCCEYDFT_Contig1-processed-gene-17.0-mRNA-1 gene=augustus-ANCCEYDFT_Contig1-processed-gene-17.0

if ( (! -e $cds2gene) or (! -e $cds_loj) ) {
    die 'Format: umass_v1.0.cds_loj_umass_v2.1.cds_genemap_06dec2022a.pl [cds2gene for *both* CDS sets] [cds_loj] > [header, then: gene1 \t gene2]', "\n";
}

my $header = "Gene\tMapped_gene";

open my $CDS2GENE, '<', $cds2gene;
while (my $input = <$CDS2GENE>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
        my $cds  = $1;
        my $gene = $2;
        if ( exists $data_ref->{'cds'}->{ $cds }->{'gene'} ) {
            die "In cds2gene table $cds2gene, redundant CDS: $cds\n";
        }
        $data_ref->{'cds'}->{ $cds }->{'gene'} = $gene;
    }
    else {
        die "From cds2gene $cds2gene, cannot parse input: $input\n";
    }
}
close $CDS2GENE;

# Format of input $cds_loj: excise two columns from a BEDtools loj join using "cut -f 10,20".
# Note heterogeneous formatting in column 2 of input.
# Sample input lines:
# 
# ID=cds:Acey_s0001.g2.t1;Parent=transcript:Acey_s0001.g2.t1      ID=cds:Acey_s0001.g2.t1;Parent=transcript:Acey_s0001.g2.t1
# 
# ID=cds:Acey_s0001.g3.t3;Parent=transcript:Acey_s0001.g3.t3      
#    ID=Acey_s0001.v2.g17938.t1.CDS1 gene=Acey_s0001.v2.g17938;Parent=Acey_s0001.v2.g17938.t1 gene=Acey_s0001.v2.g17938;

open my $CDS_LOJ, '<', $cds_loj;
while (my $input = <$CDS_LOJ>) {
    chomp $input;
    if ( $input =~ /\A [^\t]* Parent [=] (?:transcript:){0,1} (\S+) [^\t]* \t [^\t]* Parent [=] (?:transcript:){0,1} (\S+) [^\t]* \z/xms ) {
        my $cds1  = $1;
        my $gene1 = q{};
        my $cds2  = $2;
        my $gene2 = q{};

        # In the next two lines of code, reverse the order of cds1 => gene2 and cds2 => gene1, 
        #     so that I can reuse all of the subsequent code without changes,
        #     and get a reverse-order gene-centric annotation table.

        $gene2 = map_cds2gene($cds1);
        $gene1 = map_cds2gene($cds2);

        if ( ( $gene1 =~ /\S/xms ) and ( $gene2 =~ /\S/xms ) ) {
            $data_ref->{'gene1'}->{ $gene1 }->{'gene2'}->{ $gene2 } = 1;
        }
    }
    elsif ( $input !~ /\A [^\t]* Parent [=] (?:transcript:){0,1} (\S+) [^\t]* \t [.] \z/xms) {
        die "From cds_loj $cds_loj, cannot parse input: $input\n";
    }
}
close $CDS_LOJ;

my @genes1 = grep { /\S/ } sort keys %{ $data_ref->{'gene1'} };
foreach my $gene1 (@genes1) {
    my @genes2 = grep { /\S/ } sort keys %{ $data_ref->{'gene1'}->{ $gene1 }->{'gene2'} };
    my $genes2_text = join '; ', @genes2;

    print "$header\n" if $header;
    $header = q{};

    print "$gene1\t$genes2_text\n";
}

sub map_cds2gene {
    my $_cds   = $_[0];
    my $_gene  = q{};
    if ( exists $data_ref->{'cds'}->{ $_cds }->{'gene'} ) {
        $_gene = $data_ref->{'cds'}->{ $_cds }->{'gene'};
    }
    else {
        # Deal with weird habit of LiftOver: it adds '_\d+' suffixes to some genes it maps.
        $_cds =~ s/\_\d+\z//;
        if ( exists $data_ref->{'cds'}->{ $_cds }->{'gene'} ) {
            $_gene = $data_ref->{'cds'}->{ $_cds }->{'gene'};
        }
       	# Deal with variable formatting of second column:
        $_cds =~ s/\.CDS\d+\z//;
        if ( exists $data_ref->{'cds'}->{ $_cds }->{'gene'} ) {
            $_gene = $data_ref->{'cds'}->{ $_cds }->{'gene'};
       	}
        else {
            die "Cannot map CDS $_cds to a gene\n";
        }
    }
    return $_gene;
}

