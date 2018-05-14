#!/usr/bin/env perl

# wbg2omim_28apr2012.pl -- Erich Schwarz <emsch@its.caltech.edu>, 4/28/2012.
# Purpose: hacked-together data to link AWC genes to diseases.  Note that I have not really rigorously mapped everything back to WS220!

use strict;
use warnings;

# Note full WS220 gene names source was: /sternlab/redivivus/data02/schwarz/archive/taygeta/science/afd2solexa_paper/ws220_xace_data/ws220_gene_names.tsv
# Also note that WS230 humprots from WormBase *not* mapped back; nor was eggnog3 rigorously mapped!

my $gene_names  = $ARGV[0];
my $wb_humorths = $ARGV[1];
my $eggnog3     = $ARGV[2];
my $ens_omim    = $ARGV[3];

my $wbgene    = q{};
my $cds       = q{};
my $cgc       = q{};
my $full_name = q{};
my $humprot   = q{};
my $omim_id   = q{};
my $omim_txt  = q{};

my $data_ref;

if ( (! $gene_names ) or (! $wb_humorths ) or (! $eggnog3 ) or (! $ens_omim ) ) { 
    die "Format: wbg2omim_28apr2012.pl ws220_gene_names.tsv c_elegans.WS230.humprot_orths.txt eggnog3_wormhum_orths.txt ensp2omim.txt > wbgenes2omim_28apr2012.txt\n";
}

open my $GENES, '<', $gene_names or die "Can't open gene name list $gene_names: $!";
while (my $input = <$GENES>) { 
    chomp $input;

    # Example inputs:
    #   "WBGene00000024"        "AC3.3" "abu-1" "Caenorhabditis elegans"
    #   "WBGene00007063"        "2L52.1"                "Caenorhabditis elegans"

    if ( $input =~ /\A \"(WBGene\d+)\" \t \"(\S+)\" \t (\S*) \t /xms ) { 
        $wbgene    = $1;
        $cds       = $2;
        $cgc       = $3;
        $full_name = $wbgene . '|' . $cds;
        if ( $cgc =~ /\A \" (\S+) \" \z/xms ) { 
            $cgc = $1;
            $full_name .= '|' . $cgc;
        }
        $data_ref->{'cds'}->{$cds}->{'wbgene'}          = $wbgene;
        $data_ref->{'wbgene'}->{$wbgene}->{'full_name'} = $full_name;
    }
    else { 
        die "Can't parse input from $gene_names: $input\n";
    }
}
close $GENES or die "Can't close filehandle to gene name list $gene_names: $!";

open my $WB_HUMS, '<', $wb_humorths or die "Can't open WB/human ortholog list $wb_humorths: $!";
while (my $input = <$WB_HUMS>) { 
    chomp $input;
    if ( $input =~ /\A (WBGene\d+) \t (ENSP\d+) \z /xms ) { 
        $wbgene  = $1;
        $humprot = $2;
        $data_ref->{'wbgene'}->{$wbgene}->{'humorth'}->{$humprot} = 1;
        $data_ref->{'humorth'}->{$humprot}->{'wbgene'}->{$wbgene} = 1;
    }
    else { 
        die "Can't parse input from $wb_humorths: $input\n";
    }
}
close $WB_HUMS or die "Can't close filehandle to WB/human ortholog list $wb_humorths: $!";

open my $EGG, '<', $eggnog3 or die "Can't open eggnog2orth file $eggnog3: $!";
while (my $input = <$EGG>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (ENSP\d+) \z /xms ) {
        $cds     = $1;
        $humprot = $2;
        if ( exists $data_ref->{'cds'}->{$cds}->{'wbgene'} ) { 
            $wbgene  = $data_ref->{'cds'}->{$cds}->{'wbgene'};
            $data_ref->{'wbgene'}->{$wbgene}->{'humorth'}->{$humprot} = 1;
            $data_ref->{'humorth'}->{$humprot}->{'wbgene'}->{$wbgene} = 1;
        }
        if (! exists $data_ref->{'cds'}->{$cds}->{'wbgene'} ) {
            warn "Can't link CDS $cds to WBGene\n";
        }
    }
    else {
        die "Can't parse input from $wb_humorths: $input\n";
    }
}
close $EGG or die "Can't close filehandle to eggnog2orth file $eggnog3: $!";

open my $OMIM, '<', $ens_omim or die "Can't open ens2omim file $ens_omim: $!";
while (my $input = <$OMIM>) {
    chomp $input;
    if ( $input =~ /\A (ENSP\d+) \s+ (\d+) \s+ (.+) \s* \z /xms ) { 
        $humprot  = $1;
        $omim_id  = $2;
        $omim_txt = $3;
        if ($omim_id) { 
            my $omim_full = 'OMIM:' . $omim_id . ' [' . $omim_txt . ']';
            if ( exists $data_ref->{'humorth'}->{$humprot}->{'wbgene'} ) { 
                $wbgene = $data_ref->{'humorth'}->{$humprot}->{'wbgene'};
                foreach my $wb_omim_orth (sort keys %{ $data_ref->{'humorth'}->{$humprot}->{'wbgene'} } ) { 
                    $data_ref->{'wbgene'}->{$wb_omim_orth}->{'disease'}->{$omim_full} = 1;
                }
            }
        }
    }
    else {
        die "Can't parse input from $wb_humorths: $input\n";
    }
}
close $OMIM or die "Can't close filehandle to ens2omim file $ens_omim: $!";

print "Gene\tOMIM_2012\n";

foreach my $wbgene1 ( sort keys %{ $data_ref->{'wbgene'} } ) { 
    if (     ( exists $data_ref->{'wbgene'}->{$wbgene1}->{'disease'}   ) 
         and ( exists $data_ref->{'wbgene'}->{$wbgene1}->{'full_name'} ) ) { 
        $full_name = $data_ref->{'wbgene'}->{$wbgene1}->{'full_name'};
        my @diseases = sort keys %{ $data_ref->{'wbgene'}->{$wbgene1}->{'disease'} };
        my $disease_text = join '; ', @diseases;
        print "$full_name\t$disease_text\n";
    }
}

