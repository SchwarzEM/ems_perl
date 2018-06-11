#!/usr/bin/perl

# edit_Nworm_orthomcl_cdhit.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/23/2007.
# Purpose: convert N-worm orthomcl.out files proteins -> genes, using CD-HIT data to correct PB2801.

use strict;
use warnings;

my $input = "";
my $key_gene = "";
my $o_g = "";
my @other_genes = ();
my %cds2gene = ();

unless ($#ARGV == 5) { 
    die "Format: edit_3worm_orthomcl.pl",
        " [wormpep, or headers]",         # $ARGV[0]
        " [brig_prot2gene]",              # $ARGV[1]
        " [pre-rempep, or headers]",      # $ARGV[2]
        " [pre-pb2801pep, or headers]",   # $ARGV[3]
        " [cd-hit_file.txt]",             # $ARGV[4]
        " [orthomcl.out]\n";              # $ARGV[5]
}

open (WORMPEP, "$ARGV[0]") 
    or die "Can't open wormpep $ARGV[0]: $!";
while ($input = <WORMPEP>) { 
    chomp $input;
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

open (BRIG_PROT2GENE, "$ARGV[1]") 
    or die "Can't open briggsae protein-to-gene table $ARGV[1]: $!";
while ($input = <BRIG_PROT2GENE>) { 
    chomp $input;
    if ($input =~ /^BP:(CBP\d+)\s+(CBG\d+)/) {
        $cds2gene{"$1(briggsae)"} = "$2(briggsae)";
    }
}
close BRIG_PROT2GENE;

# Sample input:
# 
# >cr01.sctg0.wum.1004.1 Contig0f-snap.94.final
# >cr01.sctg0.wum.1004.2 Contig0f.Fgenesh_Celegans.63.final

open (REMPEP, "$ARGV[2]") 
    or die "Can't open preliminary rempep $ARGV[2]: $!";
while ($input = <REMPEP>) {
    chomp $input;
    if ($input =~ /^>((\S+)\.\d+)\s+/) {
        $cds2gene{"$1(remanei)"} = "$2(remanei)";
    }
}
close REMPEP;

open (PB2801PEP, "$ARGV[3]") 
    or die "Can't open preliminary PB2801pep $ARGV[3]: $!";
while ($input = <PB2801PEP>) {
    chomp $input;
    if ($input =~ /^>(Contig\d+\.\d+)\s*/) { 
        my $pb2801_cds = $1 . "(pb2801)";
        $cds2gene{$pb2801_cds} = "$pb2801_cds";
    }
}
close PB2801PEP;

open (CD_HIT, "$ARGV[4]") 
    or die "Can't open CD-HIT clustering of PB2801 $ARGV[4]: $!";
while ($input = <CD_HIT>) { 
    chomp $input;

# Sample input lines:
# 
# >Cluster 4
# 0       5291aa, >Contig442.013... *
# 1       116aa, >Contig3632.001... at 1:116:120:235/100%

    if ($input =~ /^>Cluster \d+/) { 
        if (! $key_gene) { 
            warn "No key gene before: $input\n";
        }
        if ($key_gene) {
            if (@other_genes) {
                foreach $o_g (@other_genes) { 
                    $cds2gene{$o_g} = "$key_gene";
                }
            }
            $cds2gene{$key_gene} = "$key_gene";  # backstop, but not strictly necessary
        }
        $key_gene = "";
        @other_genes = ();
    }
    elsif ($input =~ /(Contig\d+\.\d+).+\*\s*$/) { 
        $key_gene = $1 . "(pb2801)";
        if (! $cds2gene{$key_gene}) { 
            die "Misparsed PB2801 name: $key_gene?\n";
        }
    }
    elsif ($input =~ /(Contig\d+\.\d+)/) { 
        $o_g = $1 . "(pb2801)";
        if (! $cds2gene{$o_g}) { 
            die "Misparsed PB2801 name: $o_g?\n";
        }
        push @other_genes, $o_g;
    }
}

foreach my $o_g (@other_genes) {    # last round, after end of <CD_HIT>
    $cds2gene{$o_g} = "$key_gene";
}
$cds2gene{$key_gene} = "$key_gene";

close CD_HIT;

# Sample input:
# 
# ORTHOMCL3711(5 genes,4 taxa):    CBP24868(briggsae) Contig48.032(pb2801) F44B9.4a(elegans) F44B9.4b(elegans) cr01.sctg49.wum.86.1(remanei)

open (ORTHO, "$ARGV[5]") or die "Can't open N-worm-species OrthoMCL output $ARGV[5]: $!";
while ($input = <ORTHO>) { 
    chomp $input;
    if ($input =~ /^(ORTHOMCL\d+)\(\d+ genes,\d+ taxa\):\s+(.+)\s*$/) { 
#                   .$1.......$1.!!                  !!    .$2.
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
        print "$ortho_grp($ortho_no genes,$taxon_no taxa):\t @orthologs\n";
    }
}
