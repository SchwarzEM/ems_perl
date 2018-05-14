#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use Statistics::Descriptive;

my $data_ref;

my %prefix2spp = (
    'Acey' => 'Ancylostoma ceylanicum',
    'Amel' => 'Apis mellifera',
    'Asuu' => 'Ascaris suum',
    'Bflo' => 'Branchiostoma floridae',
    'Bmal' => 'Brugia malayi',
    'Bxyl' => 'B. xylophilus',
    'Cbri' => 'Caenorhabditis briggsae',
    'Cele' => 'Caenorhabditis elegans',
    'Cfam' => 'Canis familiaris',
    'Cgig' => 'Crassostrea gigas',
    'Cint' => 'Ciona intestinalis',
    'Ctel' => 'Capitella teleta',
    'Dimm' => 'Dirofilaria immitis',
    'Dmel' => 'Drosophila melanogaster',
    'Dnov' => 'Dasypus novemcinctus',
    'Drer' => 'Danio rerio',
    'Eeur' => 'Erinaceus europaeus',
    'Ggal' => 'Gallus gallus',
    'Hbac' => 'Heterorhabditis bacteriophora',
    'Hcon' => 'Haemonchus contortus',
    'Hpol' => 'Heligmosomoides (polygyrus) bakeri',
    'Hrob' => 'Helobdella robusta',
    'Hsap' => 'Homo sapiens',
    'Lafr' => 'Loxodonta africana',
    'Lcha' => 'Latimeria chalumnae',
    'Lgig' => 'Lottia gigantea',
    'Lloa' => 'Loa loa',
    'Lone' => 'Laxus oneistus',
    'Maur' => 'Mesocricetus auratus',
    'Mdom' => 'Monodelphis domestica',
    'Mhap' => 'Meloidogyne hapla',
    'Mmus' => 'Mus musculus',
    'Name' => 'Necator americanus',
    'Nbra' => 'Nippostrongylus brasiliensis',
    'Oana' => 'Ornithorhynchus anatinus',
    'Phum' => 'Pediculus humanus',
    'Pmar' => 'Petromyzon marinus',
    'Ppac' => 'Pristionchus pacificus',
    'Pred' => 'Panagrellus redivivus',
    'Skow' => 'Saccoglossus kowalevskii',
    'Spur' => 'Strongylocentrotus purpuratus',
    'Tcan' => 'Toxocara canis',
    'Tgut' => 'Taeniopygia guttata',
    'Tspi' => 'Trichinella spiralis',
    'Wban' => 'Wuchereria bancrofti',
    'Xtro' => 'Xenopus tropicalis',
);

while (my $input = <>) { 
    if ( $input !~ /\A [#] /xms ) {
        if ( $input =~ /\A (([A-Z][a-z]+) _\S+) (?: \s+ \S+){16} \s+ (\d+) \s* /xms ) { 
            my $protein        = $1;
            my $sp_prefix      = $2; 
            my $domain_count   = $3;
            if ( ( $domain_count != int($domain_count ) ) or ( $domain_count < 0 ) ) {
                die "Can't parse domain count ($domain_count) in: $input\n";
            }
            if ( $domain_count >= 1 ) { 
                push @{ $data_ref->{'sp_pref'}->{$sp_prefix}->{'domains_per_prot'} }, $domain_count;
                $data_ref->{'sp_pref'}->{$sp_prefix}->{'domain_count'}->{$domain_count}->{'protein_w_count'}->{$protein} = 1;
            }
        }
        else { 
             die "Can't parse input: $input\n";
        }
    }
}

my $header = "Prefix\tSpecies\tProts_w_doms\tMean\tMedian\tMin.\tMax.\tProts_w_max";

my @spp_prefixes = sort keys %{ $data_ref->{'sp_pref'} };

foreach my $sp_prefix (@spp_prefixes) { 
    if (! exists $prefix2spp{$sp_prefix} ) { 
        die "Can't parse species prefix: $sp_prefix\n";
    }
    my $full_name = $prefix2spp{$sp_prefix};
    my @domain_counts = @{ $data_ref->{'sp_pref'}->{$sp_prefix}->{'domains_per_prot'} };
    my $stat_dom_counts = Statistics::Descriptive::Full->new();
    $stat_dom_counts->add_data(@domain_counts); 
    my $prots_w_doms     = $stat_dom_counts->count();
    my $dom_count_mean   = $stat_dom_counts->mean();
    my $dom_count_median = $stat_dom_counts->median();
    my $dom_count_min    = $stat_dom_counts->min();
    my $dom_count_max    = $stat_dom_counts->max();

    # Make this more readable by cutting the sig. digits *way* down.
    $dom_count_mean = sprintf("%.2f", $dom_count_mean);

    # Get a list of every protein that had the maximum domain count.
    my @prots_w_max = sort keys %{ $data_ref->{'sp_pref'}->{$sp_prefix}->{'domain_count'}->{$dom_count_max}->{'protein_w_count'} };
    my $max_prot_count = @prots_w_max;
    if ( $max_prot_count > 5 ) {
        my $remainder_count = ($max_prot_count - 5); 
        @prots_w_max = @prots_w_max[0..4];
        push @prots_w_max, "and $remainder_count other(s)";
    }
    my $max_prot_text = join ', ', @prots_w_max;

    print "$header\n" if $header;
    $header = q{};

    print "$sp_prefix\t$full_name\t$prots_w_doms\t$dom_count_mean\t$dom_count_median\t$dom_count_min\t$dom_count_max\t$max_prot_text\n";
}

