#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

my %prefix2class = (
    Acey_ => 'Nematodes',
    Acpis_ => 'Arthropods',
    Acyrt_ => 'Arthropods',
    Aeaeg_ => 'Arthropods',
    Amel_ => 'Arthropods',
    Angam_ => 'Arthropods',
    Apis_ => 'Arthropods',
    Asuu_ => 'Nematodes',
    Atcep_ => 'Arthropods',
    Bflo_ => 'Metazoa',
    Bmal_ => 'Nematodes',
    Bomor_ => 'Arthropods',
    Bxyl_ => 'Nematodes',
    Cbri_ => 'Nematodes',
    Cele_ => 'Nematodes',
    Cerat_ => 'Arthropods',
    Cfam_ => 'Vertebrates',
    Cgig_ => 'Metazoa',
    Cint_ => 'Metazoa',
    Ctel_ => 'Metazoa',
    Cuqui_ => 'Arthropods',
    Daple_ => 'Arthropods',
    Dapul_ => 'Arthropods',
    Depon_ => 'Arthropods',
    Diaph_ => 'Arthropods',
    Dimm_ => 'Nematodes',
    Dmel_ => 'Arthropods',
    Dnov_ => 'Vertebrates',
    Drer_ => 'Vertebrates',
    Eeur_ => 'Vertebrates',
    Ggal_ => 'Vertebrates',
    Hbac_ => 'Nematodes',
    Hcon_ => 'Nematodes',
    Hemel_ => 'Arthropods',
    Hpol_ => 'Nematodes',
    Hrob_ => 'Metazoa',
    Hsap_ => 'Vertebrates',
    Hydra_ => 'Metazoa',
    Ixsca_ => 'Arthropods',
    Lafr_ => 'Vertebrates',
    Lcha_ => 'Vertebrates',
    Lgig_ => 'Metazoa',
    Lloa_ => 'Nematodes',
    Lone_ => 'Nematodes',
    Maur_ => 'Vertebrates',
    Mdom_ => 'Vertebrates',
    Mesca_ => 'Arthropods',
    Mhap_ => 'Nematodes',
    Mmus_ => 'Vertebrates',
    Name_ => 'Nematodes',
    Navit_ => 'Arthropods',
    Nbra_ => 'Nematodes',
    Oana_ => 'Vertebrates',
    Phum_ => 'Arthropods',
    Plano_ => 'Arthropods',
    Pmar_ => 'Vertebrates',
    Ppac_ => 'Nematodes',
    Pred_ => 'Nematodes',
    Rhpro_ => 'Arthropods',
    Skow_ => 'Metazoa',
    Soinv_ => 'Arthropods',
    Spur_ => 'Metazoa',
    Stmar_ => 'Arthropods',
    Tcan_ => 'Nematodes',
    Teurt_ => 'Arthropods',
    Tgut_ => 'Vertebrates',
    Trcas_ => 'Arthropods',
    Trich_ => 'Metazoa',
    Tspi_ => 'Nematodes',
    Wban_ => 'Nematodes',
    Xtro_ => 'Vertebrates',
);

my %prefix2secd_prefix = (
    Acey_ => 'Nema_',
    Acpis_ => 'Arth_',
    Acyrt_ => 'Arth_',
    Aeaeg_ => 'Arth_',
    Amel_ => 'Arth_',
    Angam_ => 'Arth_',
    Apis_ => 'Arth_',
    Asuu_ => 'Nema_',
    Atcep_ => 'Arth_',
    Bflo_ => 'Mzoa_',
    Bmal_ => 'Nema_',
    Bomor_ => 'Arth_',
    Bxyl_ => 'Nema_',
    Cbri_ => 'Nema_',
    Cele_ => 'Nema_',
    Cerat_ => 'Arth_',
    Cfam_ => 'Vert_',
    Cgig_ => 'Mzoa_',
    Cint_ => 'Mzoa_',
    Ctel_ => 'Mzoa_',
    Cuqui_ => 'Arth_',
    Daple_ => 'Arth_',
    Dapul_ => 'Arth_',
    Depon_ => 'Arth_',
    Diaph_ => 'Arth_',
    Dimm_ => 'Nema_',
    Dmel_ => 'Arth_',
    Dnov_ => 'Vert_',
    Drer_ => 'Vert_',
    Eeur_ => 'Vert_',
    Ggal_ => 'Vert_',
    Hbac_ => 'Nema_',
    Hcon_ => 'Nema_',
    Hemel_ => 'Arth_',
    Hpol_ => 'Nema_',
    Hrob_ => 'Mzoa_',
    Hsap_ => 'Vert_',
    Hydra_ => 'Mzoa_',
    Ixsca_ => 'Arth_',
    Lafr_ => 'Vert_',
    Lcha_ => 'Vert_',
    Lgig_ => 'Mzoa_',
    Lloa_ => 'Nema_',
    Lone_ => 'Nema_',
    Maur_ => 'Vert_',
    Mdom_ => 'Vert_',
    Mesca_ => 'Arth_',
    Mhap_ => 'Nema_',
    Mmus_ => 'Vert_',
    Name_ => 'Nema_',
    Navit_ => 'Arth_',
    Nbra_ => 'Nema_',
    Oana_ => 'Vert_',
    Phum_ => 'Arth_',
    Plano_ => 'Arth_',
    Pmar_ => 'Vert_',
    Ppac_ => 'Nema_',
    Pred_ => 'Nema_',
    Rhpro_ => 'Arth_',
    Skow_ => 'Mzoa_',
    Soinv_ => 'Arth_',
    Spur_ => 'Mzoa_',
    Stmar_ => 'Arth_',
    Tcan_ => 'Nema_',
    Teurt_ => 'Arth_',
    Tgut_ => 'Vert_',
    Trcas_ => 'Arth_',
    Trich_ => 'Mzoa_',
    Tspi_ => 'Nema_',
    Wban_ => 'Nema_',
    Xtro_ => 'Vert_',
);

my %prefix2species = (
    Acey_ => 'Ancylostoma ceylanicum',
    Amel_ => 'Apis mellifera',
    Asuu_ => 'Ascaris suum',
    Bflo_ => 'Branchiostoma floridae',
    Bmal_ => 'Brugia malayi',
    Bxyl_ => 'Bursaphelenchus xylophilus',
    Cbri_ => 'Caenorhabditis briggsae',
    Cele_ => 'Caenorhabditis elegans',
    Cfam_ => 'Canis familiaris',
    Cgig_ => 'Crassostrea_gigas',
    Cint_ => 'Ciona intestinalis',
    Ctel_ => 'Capitella teleta',
    Dimm_ => 'Dirofilaria immitis',
    Dmel_ => 'Drosophila melanogaster',
    Dnov_ => 'Dasypus novemcinctus',
    Drer_ => 'Danio rerio',
    Eeur_ => 'Erinaceus europaeus',
    Ggal_ => 'Gallus gallus',
    Hbac_ => 'Heterorhabditis bacteriophora',
    Hcon_ => 'Haemonchus contortus',
    Hpol_ => 'Heligmosomoides polygyrus',
    Hrob_ => 'Helobdella robusta',
    Hsap_ => 'Homo sapiens',
    Lafr_ => 'Loxodonta africana',
    Lcha_ => 'Latimeria chalumnae',
    Lgig_ => 'Lottia gigantea',
    Lloa_ => 'Loa loa',
    Lone_ => 'Laxus oneistus',
    Maur_ => 'Mesocricetus auratus',
    Mdom_ => 'Monodelphis domestica',
    Mhap_ => 'Meloidogyne hapla',
    Mmus_ => 'Mus musculus',
    Name_ => 'Necator americanus',
    Nbra_ => 'Nippostrongylus brasiliensis',
    Oana_ => 'Ornithorhynchus anatinus',
    Phum_ => 'Pediculus humanus',
    Pmar_ => 'Petromyzon marinus',
    Ppac_ => 'Pristionchus pacificus',
    Pred_ => 'Panagrellus redivivus',
    Skow_ => 'Saccoglossus kowalevskii',
    Spur_ => 'Strongylocentrotus purpuratus',
    Tcan_ => 'Toxocara canis',
    Tgut_ => 'Taeniopygia guttata',
    Tspi_ => 'Trichinella spiralis',
    Wban_ => 'Wuchereria bancrofti',
    Xtro_ => 'Xenopus tropicalis',
);

while (my $input = <> ) {
    chomp $input;
    if ( $input =~ /\A > /xms ) {
        my $seqname = q{};
        my $prefix  = q{};
        my $protein = q{};

        if ( $input =~ /\A > (([A-Za-z]+ [_]) (\S+)) /xms ) { 
            $seqname = $1;
            $prefix  = $2;
            $protein = $3;
            if (    (! exists $prefix2class{$prefix}       ) 
                 or (! exists $prefix2species{$prefix}     ) 
                 or (! exists $prefix2secd_prefix{$prefix} ) ) {
                die "Cannot parse prefix \"$prefix\" in: $input\n";
            }
            else {
                my $class    = $prefix2class{$prefix};
                my $species  = $prefix2species{$prefix};
                my $sec_pref = $prefix2secd_prefix{$prefix};
                my $gene     = q{};

# Sample ENSEMBL headers:
#
# >Amel_GB10001  GB10001-PA pep:novel supercontig:Amel4.0:Group2.31:545851:547406:1 gene:GB10001
# >Cfam_ENSCAFG00000000001  ENSCAFP00000000001 pep:known_by_projection chromosome:CanFam3.1:1:247829:322180:-1 gene:ENSCAFG00000000001
# >Ctel_CapteG100001  CapteP100001 pep:novel supercontig:GCA_000328365.1:CAPTEscaffold_335:206104:207618:-1 gene:CapteG100001
# >Cint_ENSCING00000000016  ENSCINP00000000022 pep:known chromosome:KH:2:1522050:1524495:-1 gene:ENSCING00000000016
# >Cgig_CGI_10000009  EKC17950 pep:novel supercontig:GCA_000297895.1:C18346:174:551:1 gene:CGI_10000009
# >Drer_ENSDARG00000000102  ENSDARP00000000102 pep:known chromosome:Zv9:10:13048113:13142710:-1 gene:ENSDARG00000000102
# >Dnov_ENSDNOG00000000007  ENSDNOP00000000004 pep:known_by_projection scaffold:Dasnov3.0:JH568938.1:1313231:1387790:1 gene:ENSDNOG00000000007
# >Dmel_FBgn0000015  FBpp0082826 pep:known chromosome:BDGP5:3R:12753395:12760298:-1 gene:FBgn0000015
# >Eeur_ENSEEUG00000000001  ENSEEUP00000000001 pep:novel scaffold:HEDGEHOG:scaffold_327592:5923:6423:1 gene:ENSEEUG00000000001
# >Ggal_ENSGALG00000000003  ENSGALP00000000002 pep:known_by_projection chromosome:Galgal4:1:19947199:19967662:1 gene:ENSGALG00000000003
# >Hrob_HelroG100002  HelroP100002 pep:novel supercontig:GCA_000326865.1:HELROscaffold_26:1630662:1637446:-1 gene:HelroG100002
# >Hsap_ENSG00000000003  ENSP00000362111 pep:known chromosome:GRCh37:X:99883667:99891803:-1 gene:ENSG00000000003
# >Lcha_ENSLACG00000000001  ENSLACP00000000001 pep:novel scaffold:LatCha1:JH131052.1:2:6373:-1 gene:ENSLACG00000000001
# >Lloa_EFO28458.1  EFO28458.1	EFO28458.1	gene:LOAG_00002
# >Lgig_LotgiG100015  LotgiP100015 pep:novel supercontig:GCA_000327385.1:LOTGIsca_37:2203307:2203474:1 gene:LotgiG100015
# >Lafr_ENSLAFG00000000001  ENSLAFP00000000001 pep:known_by_projection supercontig:loxAfr3:scaffold_36:2062840:2100909:-1 gene:ENSLAFG00000000001
# >Mdom_ENSMODG00000000001  ENSMODP00000000002 pep:novel chromosome:BROADO5:7:122757216:122916959:-1 gene:ENSMODG00000000001
# >Mmus_ENSMUSG00000000001  ENSMUSP00000000001 pep:known chromosome:GRCm38:3:108107280:108146146:-1 gene:ENSMUSG00000000001
# >Oana_ENSOANG00000000110  ENSOANP00000000222 pep:known_by_projection supercontig:OANA5:Contig20092:18594:22750:1 gene:ENSOANG00000000110
# >Phum_PHUM000210  PHUM000210-PA pep:known supercontig:PhumU1:DS234986:20136:32801:1 gene:PHUM000210
# >Pmar_ENSPMAG00000000001  ENSPMAP00000000001 pep:known_by_projection scaffold:Pmarinus_7.0:GL476426:6020:9121:-1 gene:ENSPMAG00000000001
# >Tgut_ENSTGUG00000000001  ENSTGUP00000000001 pep:known_by_projection chromosome:taeGut3.2.4:24:16660:23292:-1 gene:ENSTGUG00000000001
# >Tspi_EFV58942  EFV58942	EFV58942	gene:Tsp_00001
# >Xtro_ENSXETG00000000002  ENSXETP00000000002 pep:known scaffold:JGI_4.2:GL173505.1:64854:70530:1 gene:ENSXETG00000000002

                # ENSEMBL-derived header line syntax:
                if ( $input =~ /\A > \S+ \s+ (\S+) \s .+ \s gene:(\S+) \s* \z/xms ) {
                    $protein = $1;
                    $gene    = $2;
                }

# Sample NCBI headers:
# 
# >Bflo_EEN41500.1  gi|229270477|gb|EEN41500.1| hypothetical protein BRAFLDRAFT_111961 [Branchiostoma floridae]
# >Tcan_AAD31000.1  gi|4838459|gb|AAD31000.1|AF126830_1 C-type lectin Tc-ctl-4 [Toxocara canis]
# >Tcan_AAB96779.1  gi|2773355|gb|AAB96779.1| excretory/secretory C-type lectin TES-32 [Toxocara canis]
# >Lone_AAX22004.1  gi|60499547|gb|AAX22004.1| mermaid-1 [Laxus oneistus]
# >Hpol_ACS37721.1  gi|240002899|gb|ACS37721.1| C-type lectin-1 [Heligmosomoides polygyrus]
# >Nbra_ACS37722.1  gi|240002901|gb|ACS37722.1| C-type lectin-1 [Nippostrongylus brasiliensis]
# >Nbra_ACS37723.1  gi|240002903|gb|ACS37723.1| C-type lectin-2 [Nippostrongylus brasiliensis]
# >Tcan_BAA31253.2  gi|4586556|dbj|BAA31253.2| proteoglycan core protein [Toxocara canis]
# >Maur_XP_005063187.1  gi|524917775|ref|XP_005063187.1| PREDICTED: GDNF family receptor alpha-1 [Mesocricetus auratus]
# >Skow_XP_006811105.1  gi|585641944|ref|XP_006811105.1| PREDICTED: laminin subunit alpha-2-like [Saccoglossus kowalevskii]
# >Spur_XP_003723254.1  gi|390331338|ref|XP_003723254.1| PREDICTED: uncharacterized protein LOC100887917 [Strongylocentrotus purpuratus]
# >Wban_EJW69766.1  gi|402575806|gb|EJW69766.1| hypothetical protein WUBG_19327 [Wuchereria bancrofti]

                # NCBI-derived header line syntax:
                elsif ( $input =~ /\A > \S+ \s+ gi\|\d+ /xms ) {
                    $gene = $protein;
                }

# Sample WBgene headers:
# 
# >Bmal_Bm14  Bm14        BM26371 WBGene00220275  status:Predicted
# >Cbri_CBG00014  CBG00014        CBP36006        WBGene00023529  status:Predicted        UniProt:A8WM50
# >Cele_abu-3  F31A3.1    CE07158 WBGene00000026  locus:abu-3     status:Predicted        UniProt:Q19919  protein_id:CCD70290.1
# >Ppac_PPA00001  PPA00001	PP37808	WBGene00089555	status:Predicted

                # WormBase-derived  header line syntax:
                elsif ( $input =~ /\A > \S+ \s+ \S+ \s+ (\S+) \s+ (WBGene\d+) /xms ) {
                    $protein = $1;
                    $gene    = $2;
                }

# Sample basic nematode headers:
# 
# >Asuu_ASU_00021  ASU_00021      ASU_00021       ASU_00021
# >Bxyl_BUX.c00008.1  BUX.c00008.1        BUX.c00008.1    BUX.gene.c00008.1
# >Dimm_nDi.2.2.2.t00163  nDi.2.2.2.t00163	nDi.2.2.2.t00163	nDi.2.2.2.g00163
# >Hbac_Hba_00001  Hba_00001	Hba_00001	Hba_00001
# >Hcon_HCOI00000100.t1  HCOI00000100.t1	HCOI00000100.t1	HCOI00000100
# >Mhap_MhA1_Contig0.frz3.gene72  MhA1_Contig0.frz3.gene72	MhA1_Contig0.frz3.gene72	MhA1_Contig0.frz3.gene72
# >Name_NECAME_11531  NECAME_11531	NECAME_11531	NECAME_11531
# >Pred_g1.t1  g1.t1	g1.t1	g1

                # Really basic syntax of some nematode proteomes:
                elsif ( $input =~ /\A > \S+ \s+ \S+ \s+ (\S+) \s+ (\S+) \s* \z/xms ) {
                    $protein = $1;
                    $gene    = $2;
                }

# Sample Acey header:
#
# >Acey_s0001.g1.t1

                # Acey:
                elsif ( $input =~ /\A > ((Acey_s\d+\.g\d+)\.t\d+) \b/xms ) {
                    $protein = $1;
                    $gene    = $2;
                }

# Oddball exception:
# >Asuu_MRCL   Asuu_ASU_08002  ASU_08002        ASU_08002       ASU_08002
                
                elsif ( $input =~ /\A > Asuu_MRCL \s+ \S+ \s+ \S+ \s+ (\S+) \s+ (\S+) \s* \z/xms ) {
                    $protein = $1;
                    $gene    = $2;
                }

# Other oddball exceptions:
# >Bmal_MRCLpart1  MRCLpart1_BMAL
# >Bmal_MRCLpart2   MRCLpart2_BMAL

                elsif ( $input =~ /\A > (Bmal_MRCLpart\d+) \s /xms ) {
                    $protein = $1;
                    $gene    = $protein;
                }

# Yet others:
# >Dimm_MRCL1  MRCL1_DIMM  [from: exonerate ...
# >Lloa_MRCL1  MRCL1_LOA ...

                elsif ( $input =~ /\A > ([A-Z][a-z]+_MRCL1) \s /xms ) {
                    $protein = $1;
                    $gene    = $protein
                }

                # Enforce unambiguous parsing:
                else {
                    die "Failed to parse protein header line: $input\n";
                }

                # Finally:
                print "$sec_pref$seqname\t$class\t$species\t$protein\t$gene\n";
            }
        }
        else { 
            die "Cannot parse sequence name in: $input\n";
        }
    }
}

