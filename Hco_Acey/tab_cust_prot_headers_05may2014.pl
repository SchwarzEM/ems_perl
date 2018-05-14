#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my $header = "Name\tGene\tSpecies\tSource\tDetails\n";

my %prefix2species = (
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
    'Hpol' => 'Heligmosomoides polygyrus',
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

my %prefix2source = (
    'Acey' => 'This study',
    'Amel' => 'Ensembl Genomes 21',
    'Asuu' => 'WormBase WS242',
    'Bflo' => 'NCBI',
    'Bmal' => 'WormBase WS242',
    'Bxyl' => 'WormBase WS242',
    'Cbri' => 'WormBase WS242',
    'Cele' => 'WormBase WS242',
    'Cfam' => 'Ensembl 75',
    'Cgig' => 'Ensembl Genomes 21',
    'Cint' => 'Ensembl 75',
    'Ctel' => 'Ensembl Genomes 21',
    'Dimm' => 'WormBase WS242',
    'Dmel' => 'Ensembl 75',
    'Dnov' => 'Ensembl 75',
    'Drer' => 'Ensembl 75',
    'Eeur' => 'Ensembl 75',
    'Ggal' => 'Ensembl 75',
    'Hbac' => 'WormBase WS242',
    'Hcon' => 'WormBase WS242',
    'Hpol' => 'NCBI',
    'Hrob' => 'Ensembl Genomes 21',
    'Hsap' => 'Ensembl 75',
    'Lafr' => 'Ensembl 75',
    'Lcha' => 'Ensembl 75',
    'Lgig' => 'Ensembl Genomes 21',
    'Lloa' => 'WormBase WS242',
    'Lone' => 'NCBI',
    'Maur' => 'NCBI',
    'Mdom' => 'Ensembl 75',
    'Mhap' => 'WormBase WS242',
    'Mmus' => 'Ensembl 75',
    'Name' => 'WormBase WS242',
    'Nbra' => 'NCBI',
    'Oana' => 'Ensembl 75',
    'Phum' => 'Ensembl Genomes 21',
    'Pmar' => 'Ensembl 75',
    'Ppac' => 'WormBase WS242',
    'Pred' => 'WormBase WS242',
    'Skow' => 'NCBI',
    'Spur' => 'NCBI',
    'Tcan' => 'NCBI',
    'Tgut' => 'Ensembl 75',
    'Tspi' => 'WormBase WS242',
    'Wban' => 'NCBI',
    'Xtro' => 'Ensembl 75', 
);

my $data_ref;

while (my $input = <>) { 
    chomp $input;
    if ( $input =~ /\A > /xms ) { 
        if ( $input =~ /\A > (([A-Za-z0-9]+) _ (\S+)) \b (.*) \z/xms ) {
            my $name   = $1;
            my $prefix = $2;
            my $gene   = $3;
            my $notes  = $4;
            if (! exists $prefix2species{$prefix} ) { 
                die "Cannot identify species for prefix $prefix\n";
            }
            if (! exists $prefix2source{$prefix} ) {
                die "Cannot identify data source for prefix $prefix\n";
            }
            if ( exists $data_ref->{'prefix'}->{$prefix}->{'name'}->{$name} ) { 
                die "Redundant sequence (\"$name\") in header: $input\n";
            }
            $notes =~ s/\A\s+//;
            $notes =~ s/\s+\z//;
            $data_ref->{'prefix'}->{$prefix}->{'name'}->{$name}->{'notes'} = $notes;
            $data_ref->{'prefix'}->{$prefix}->{'name'}->{$name}->{'gene'} = $gene;
        }
        else { 
            die "Can't parse input header line: $input\n";
        }
    }
}

my @prefixes = sort keys %{ $data_ref->{'prefix'} };
foreach my $prefix (@prefixes) {
    my $species = $prefix2species{$prefix};
    my $source  = $prefix2source{$prefix};
    my @names = sort keys %{ $data_ref->{'prefix'}->{$prefix}->{'name'} };
    foreach my $name (@names) { 
        my $notes = q{};
        if ( exists $data_ref->{'prefix'}->{$prefix}->{'name'}->{$name}->{'notes'} ) {
            $notes = $data_ref->{'prefix'}->{$prefix}->{'name'}->{$name}->{'notes'};
        }
        my $gene = $data_ref->{'prefix'}->{$prefix}->{'name'}->{$name}->{'gene'};
        # This should only happen once, at the top of the output, *if* there is anything to output:
        print $header if $header;
        $header = q{};

        print "$name\t$gene\t$species\t$source\t$notes\n";
    }
}

