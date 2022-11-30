#!/usr/bin/env perl

use strict ;
use warnings;
use autodie;

my $cds2gene = q{};
my $cds_loj  = q{};

$cds2gene = $ARGV[0] if $ARGV[0];
$cds_loj  = $ARGV[1] if $ARGV[1];

my $header = "Gene\tMapped_gene";

my $data_ref;   # declare as completely unspecified variable, *not* as ' = q{}', or data structures will crash under 'use strict'.

# Format of input $cds2gene: "CDS/transcript name \t gene name"
# Extract this from, e.g., ParaSite protein FASTA header lines such as:
# 
# >augustus-ANCCEYDFT_Contig1-processed-gene-17.0-mRNA-1 transcript=augustus-ANCCEYDFT_Contig1-processed-gene-17.0-mRNA-1 gene=augustus-ANCCEYDFT_Contig1-processed-gene-17.0

if ( (! -e $cds2gene) or (! -e $cds_loj) ) {
    die 'Format: washu.cds_loj_umass.cds_genemap_30nov2022b.pl [cds2gene for *both* CDS sets] [cds_loj] > [annot table: gene2 \t gene1]', "\n";
}

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
# ID=cds:maker-ANCCEYDFT_Contig2418-pred_gff_fgenesh-gene-0.0-mRNA-1;Parent=transcript:maker-ANCCEYDFT_Contig2418-pred_gff_fgenesh-gene-0.0-mRNA-1;extra_copy_number=0	
#     ID=Acey_s1646.v2.g4956.t1.CDS3 gene=Acey_s1646.v2.g4956;Parent=Acey_s1646.v2.g4956.t1 gene=Acey_s1646.v2.g4956;
# 
# ID=cds:maker-ANCCEYDFT_Contig1171-pred_gff_snap-gene-0.3-mRNA-1;Parent=transcript:maker-ANCCEYDFT_Contig1171-pred_gff_snap-gene-0.3-mRNA-1;extra_copy_number=0	
#     ID=cds:Acey_s1656.g3934.t1;Parent=transcript:Acey_s1656.g3934.t1
# 
# ID=cds:maker-ANCCEYDFT_Contig1171-pred_gff_snap-mRNA-0.4-mRNA-1;Parent=transcript:maker-ANCCEYDFT_Contig1171-pred_gff_snap-mRNA-0.4-mRNA-1;extra_copy_number=0	
#     ID=Acey_s1705.v2.g2.t1.CDS2 gene=Acey_s1705.v2.g2;Parent=Acey_s1705.v2.g2.t1 gene=Acey_s1705.v2.g2;

open my $CDS_LOJ, '<', $cds_loj;
while (my $input = <$CDS_LOJ>) {
    chomp $input;
    if ( $input =~ /\A ID[=]cds[:] ([^\s;]+) .* \t ID[=](?:cds[:]){0,1} ([^\s;]+) .* \z/xms ) {
        my $cds1  = $1;
        my $gene1 = q{};
        my $cds2  = $2;
        my $gene2 = q{};

        # In these two lines, link $cds1 to $gene2, not $gene1 -- and vice versa.
        # This lets the remainder of the code be reused with almost no changes.
        $gene2 = map_cds2gene($cds1);
        $gene1 = map_cds2gene($cds2);

        if ( ( $gene1 =~ /\S/xms ) and ( $gene2 =~ /\S/xms ) ) {
            $data_ref->{'gene1'}->{ $gene1 }->{'gene2'}->{ $gene2 } = 1;
        }
    }
    elsif ( $input !~ /\A .+ \t [.] \z/xms) {
        die "From cds_loj $cds_loj, cannot parse input: $input\n";
    }
}
close $CDS_LOJ;

my @genes1 = grep { /\S/ } sort keys %{ $data_ref->{'gene1'} };
foreach my $gene1 (@genes1) {
    my @genes2 = grep { /\S/ } sort keys %{ $data_ref->{'gene1'}->{ $gene1 }->{'gene2'} };
    my $gene2_text = join '; ', @genes2;
    print "$header\n" if $header;
    $header = q{};
    print "$gene1\t$gene2_text\n";
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

