#!/usr/bin/perl

# aceify_HymanRNAis.pl: Erich Schwarz <emsch@its.caltech.edu>, 5/3/05.
# Purpose: de novo .ace-ification of Hyman's RNAi table.  (Does not wipe older phenotypes!)

unless ($#ARGV == 1) { die "Format:  aceify_HymanRNAis.pl  [Hyman's table]  [Tablemaker RNAi list]\n"; }

# The expected (tab-delimited) format for 'Hyman's table' is:
# 
# "B0365.7"	"502505"	"2005390"	"137-b2"	
#     "TAATACGACTCACTATAGGAGACTACAGCCGCAACTGGT"	"AATTAACCCTCACTAAAGGCAAGCCTTTCAGCTCTTTCC"	
#     "Osmotic Integrity"	"Osmotic sensitivity observed in 2 embryos, mostly wildtype 
#      recordings for the remainder."	"Wild type"	"3/3 tests are wildtype"

# The expected (tab-delimited) format for 'Tablemaker RNAi list' is:
# 
# "WBRNAi00038987"	"TH:137-b2"	"WT"
#      RNAi          Historical name  Phenotype

$hyman_table = $ARGV[0];
$rnai_list   = $ARGV[1];

open (HYMAN, "$hyman_table");  # First, get circumscribed list of RNAis to ace-ify.
while (<HYMAN>) { 
    chomp ($input = $_);
    if ($input =~ /^(\"[^\"]+\"\s+){3}\"([^\"]+)\"\s+/) { 
        $hyman_id = "TH:" . $2;
        $hyman_names{$hyman_id} = 1;
    }
}
close HYMAN;

open (WBASE, "$rnai_list");    # Second, store correct names for RNAis.

while (<WBASE>) { 
    chomp ($input = $_);
    if ($input =~ /^\"([^\"]+)\"\s+\"([^\"]+)\"\s+\"([^\"]+)\"/) { 
        $wbase_rna_id = $1;
        $hyman_id     = $2;
        $old_pheno    = $3;
        if ($hyman_names{$hyman_id}) { 
            $hyman_names{$hyman_id} = $wbase_rna_id;
        }
    }
}
close WBASE;

%correct_phenos = ( 

# Hyman's name:                              # Corresponding WormBase term:

'Aberrant Cytoplasmic Structures',           'Aberrant Cytoplasmic Structures',
'Adult Phenotype-Dumpy (dpy)',               'Dpy',
'Adult Phenotype-Egg Laying Defect (egl)',   'Egl',
'Adult Phenotype-Protruding Vulva (pvl)',    'Pvl',
'Adult Phenotype-Roller (rol)',              'Rol',
'Adult Phenotype-Uncoordinated (unc)',       'Unc',
'Asymmetry of Division',                     'Asymmetry of Division abnormal',
'Centrosome Attachment',                     'Centrosome Attachment abnormal',
'Chromosome Segregation (karyomeres)',       'Chromosome Segregation (karyomeres) abnormal',
'Complex Phenotype',                         'Complex Phenotype',
'Cortical Dynamics',                         'Cortical Dynamics abnormal',
'Cytokinesis',                               'Cyk',
'Developmental Delay',                       'Developmental Delay',
'Egg Size',                                  'Egg Size abnormal',
'Embryonic Lethal',                          'Emb',
'Entry Into Interphase',                     'Entry Into Interphase abnormal',
'General Pace of Development',               'General Pace of Development abnormal',
'Integrity of Membranous Organelles',        'Integrity of Membranous Organelles defective',
'Larval Arrest-Early (L1/L2)',               'Larval Arrest-Early (L1/L2)',
'Larval Arrest-Late (L3/L4)',                'Larval Arrest-Late (L3/L4)',
'Larval Lethal-Early (L1/L2)',               'Larval Lethal-Early (L1/L2)',
'Larval Lethal-Late (L3/L4)',                'Larval Lethal-Late (L3/L4)',
'Larval-Specific Phenotype',                 'Larval-Specific Phenotype',
'Nuclear Appearance',                        'Nuclear Appearance abnormal',
'Osmotic Integrity',                         'Osmotic Integrity defective',
'Pace of P-Lineage',                         'Pace of P-Lineage abnormal',
'Passage Through Meiosis',                   'Passage Through Meiosis defective',
'Polar Body Extrusion',                      'Passage Through Meiosis defective',
'Pronuclear Migration',                      'Pnm',
'Pronuclear/Nuclear Appearance',             'Pronuclear/Nuclear Appearance abnormal',
'Severe Pleiotropic Defects',                'Severe Pleiotropic Defects',
'Sister Chromatid Separation (Cross-eyed)',  'Sister Chromatid Separation abnormal (Cross-eyed)',
'Spindle Assembly',                          'Spindle Assembly abnormal',
'Spindle Elongation/Integrity',              'Spindle Elongation/Integrity abnormal',
'Sterile F0/Fertility Problems',             'Sterile F0/Fertility Problems',
'Sterile F1',                                'Sterile F1',
'Sterility/Impaired Fertility in F0',        'Sterility/Impaired Fertility in F0', 
'Wild type',                                 'WT', 
);

open (HYMAN, "$hyman_table");  # Finally, crank out the .ace.
print "\n";
print "// NEW annotations for Hyman RNAi phenotypes.\n";
print "\n";

# sub rnai_prescoreds {  # check if $correct_phenos{$emb_phen} is element of @{$prescored{$hyman_id}} 
#     my ($query, $prescored) = @_;
#     foreach $p (@{$prescored}) { 
#         if ($query == $p) { 
#             return 0;  # kinda buggy
#         }
#     }
# return 1;
# }

while (<HYMAN>) { 
    chomp ($input = $_);
    $rnai = ""; $emb_phen = ""; $postemb_phen = ""; 
    if ($input =~ /^(\"[^\"]*?\"\s+){3}\"([^\"]*?)\"\s+(\"[^\"]*?\"\s+){2}\"([^\"]*?)\"\s+\"[^\"]*?\"\s+\"([^\"]*?)\"\s+/) { 
        $hyman_id     = "TH:" . $2;
        $emb_phen     = $4;
        $postemb_phen = $5; 
        if ($hyman_names{$hyman_id} =~ /^WBRNAi\d{8}/) { 
            if ( $correct_phenos{$emb_phen}) { 
                    push ( @{$scored_rnai{$hyman_id}}, $correct_phenos{$emb_phen}); 
                } 
                if ($correct_phenos{$postemb_phen}) { 
                    push ( @{$scored_rnai{$hyman_id}}, $correct_phenos{$postemb_phen});
                }
        } 
    } 
} 
close HYMAN;

foreach $wb_rnai (sort keys %scored_rnai) { 
    %already_scored = "";
    print "RNAi : \"$hyman_names{$wb_rnai}\"    \/\/ $wb_rnai\n";
    foreach $rnai_pheno (sort @{$scored_rnai{$wb_rnai}}) { 
        unless ($already_scored{$rnai_pheno}) { 
            print "Phenotype \"$rnai_pheno\"\n";
            $already_scored{$rnai_pheno} = 1;
        }
    }
    print "\n";
}
