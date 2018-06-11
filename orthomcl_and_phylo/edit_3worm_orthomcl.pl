#!/usr/bin/perl

# edit_3worm_orthomcl.pl -- Erich Schwarz <emsch@its.caltech.edu>, 9/10/2006.
# Purpose: convert misleading three-worm protein-based orthomcl.out files to non-misleading gene-based ones.

my %cds2gene = ();

unless ($#ARGV == 3) { die "Format: edit_3worm_orthomcl.pl [wormpep] [brig_prot2gene] [pre-rempep] [orthomcl.out]\n"; }

open (WORMPEP, "$ARGV[0]") or die "Can't open wormpep $ARGV[0]: $!";
while (<WORMPEP>) { 
    chomp($input = $_);
    if ($input =~ /^>(\S+)\s+.*\s+(WBGene\d+)\s+/) { 
        $cds2gene{"$1(elegans)"} = "$2(elegans)";
    }
}
close WORMPEP;

# Sample input:
# 
# BP:CBP00003     CBG00001
# BP:CBP20727     CBG00002
# 
# This is generated by ws_brig_prot2gene.pl from a local ACeDB database.

open (BRIG_PROT2GENE, "$ARGV[1]") or die "Can't open briggsae protein-to-gene table $ARGV[1]: $!";
while (<BRIG_PROT2GENE>) { 
    chomp($input = $_);
    if ($input =~ /^BP:(CBP\d+)\s+(CBG\d+)/) {
        $cds2gene{"$1(briggsae)"} = "$2(briggsae)";
    }
}
close BRIG_PROT2GENE;

# Sample input:
# 
# >cr01.sctg0.wum.1004.1 Contig0f-snap.94.final
# >cr01.sctg0.wum.1004.2 Contig0f.Fgenesh_Celegans.63.final

open (REMPEP, "$ARGV[2]") or die "Can't open preliminary rempep $ARGV[2]: $!";
while (<REMPEP>) {
    chomp($input = $_);
    if ($input =~ /^>((\S+)\.\d+)\s+/) {
        $cds2gene{"$1(remanei)"} = "$2(remanei)";
    }
}
close REMPEP;

# Sample input:
# 
# ORTHOMCL4096(4 genes,3 taxa):    B0336.11a(elegans) CBP02197(briggsae) cr01.sctg5.wum.325.1(remanei) cr01.sctg5.wum.325.2(remanei)
# ORTHOMCL4131(4 genes,3 taxa):    B0035.18(elegans) CBP11420(briggsae) Y111B2A.2(elegans) cr01.sctg123.wum.17.1(remanei)

open (ORTHO, "$ARGV[3]") or die "Can't open three-worm-species OrthoMCL output $ARGV[3]: $!";
while (<ORTHO>) { 
    chomp($input = $_);
    if ($input =~ /^(ORTHOMCL\d+)\(\d+ genes,\d+ taxa\):\s+(.+)\s*$/) { 
        my $ortho_grp = $1;
        my @orthoprots = split /\s+/, $2;
        my %species_seen = ();
        my %genes_seen = ();
        foreach my $o_prot (@orthoprots) { 
            if ($o_prot =~ /.+\((\w+)\)/) { $species_seen{$1} = 1; }
            $genes_seen{$cds2gene{$o_prot}} = 1;
        }
        my $ortho_no  = scalar(keys %genes_seen);
        my $taxon_no  = scalar(keys %species_seen);
        my @orthologs = sort keys %genes_seen;
        print "$ortho_grp($ortho_no genes,$taxon_no taxa)\t @orthologs\n";
    }
}
