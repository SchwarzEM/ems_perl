#!/usr/bin/perl

# wipe_past_HymanRNAi_phenos.pl
# Erich Schwarz <emsch@its.caltech.edu>, 5/3/05.
# Purpose: scan a previous .ace file of RNAis and delete *only* phenotypes for Hyman RNAis.

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

open (HYMAN, "$hyman_table");  # First, get circumscribed list of only Hyman RNAis to edit.
while (<HYMAN>) { 
    chomp ($input = $_);
    if ($input =~ /^(\"[^\"]+\"\s+){3}\"([^\"]+)\"\s+/) { 
        $hyman_id = "TH:" . $2;
        $hyman_names{$hyman_id} = 1;
    }
}
close HYMAN;

open (WBASE, "$rnai_list");    # Second, store correct names for RNAis and wipe their old phenos.
print "\n";
print "// Deleting previous Hyman RNAi phenotype annotations.\n";
print "\n";

while (<WBASE>) { 
    chomp ($input = $_);
    if ($input =~ /^\"([^\"]+)\"\s+\"([^\"]+)\"\s+\"([^\"]+)\"/) { 
        $wbase_rna_id = $1;
        $hyman_id     = $2;
        $old_pheno    = $3;
        if ($hyman_names{$hyman_id}) { 
            $hyman_names{$hyman_id} = $wbase_rna_id;
            print "RNAi : \"$wbase_rna_id\"\n";
            print "-D Phenotype \"$old_pheno\"\n\n";
        }
    }
}
close WBASE;

