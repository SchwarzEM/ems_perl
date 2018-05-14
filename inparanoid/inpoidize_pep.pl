#!/usr/bin/perl

# inpoidize_wpep.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/13/2007
# Purpose: convert (worm|brig|rem)pep headers into InParanoid-preferred headers.

use strict;
use warnings;

my $input   = "";    # input line
my $cds     = "";    # e.g., "T13A10.10a"
my $prot    = "";    # e.g., "CE30174" 
my $gene    = "";    # e.g., "WBGene00000005"
my $mgen    = "";    # e.g., "T13A10.10" in "T13A10.10a" (or "C07H6.7" in "C07H6.7")
my $cgc     = "";    # e.g., "aat-4" in "locus:aat-4"
my $blurb   = "";    # e.g., "amino acid permease"
my $uniprot = "";    # e.g., "TR:Q22448"
my $prot_id = "";    # e.g., "AAB38125.2"

unless ($#ARGV == 3) { 
    die "Format: ./inpoidize_wpep.pl [fasta] [source DB] [version] [species]\n";
}
my $species = pop @ARGV;
my $version = pop @ARGV;
my $source  = pop @ARGV;

my %ok_spp = ( elegans  => 1,
               briggsae => 1,
               remanei  => 1, );
if (! $ok_spp{$species} ) {
    die "Species name \"$species\" not recognized\n";
}

while ($input = <>) { 
    # If a non-empty non-header line, just pass on unchanged:
    if ( ($input !~ /^>/) and ($input =~ /\S/) ) { 
        print $input;
    }

    # Importing header lines from brigpep2gid (see below [1] for details):
    elsif ( ($species eq 'briggsae')
             and ($input =~ /\A > (CBP\d+)               # $1    == $prot
                             \s+ (WBGene\d+)?            # $2    == $gene
                             \|  (CBG\d+)                # $3    == $cgc
                             ( \s+ (SW|TR) : (\S+)       # $6    == $uniprot
                             )?                          # [$4]  == [optional element]
                             ( \s+ protein_id: (\S+)     # $8    == $prot_id
                             )?                          # [$7] == [optional element]
                             /x) ) {                     # ignore whitespace in regex
        $prot    = $1;
        $gene    = $2;
        $cgc     = $3;
        $uniprot = $6;
        $prot_id = $8;

        # Designed to deal with seven freak cases in brigpep2(gid) that have 
        #    a CBG\n+ name but not a WBGene\n+ name:

        unless ($gene) { 
            $gene = $cgc;
        }

        # Exporting InParanoid-ized header lines:
        print ">$prot";
        print " source=$source;";
        print " version=$version;";
        print " symbol=$cgc;"       if ($cgc);
        print " uniprot=$uniprot;"  if ($uniprot);

        print " desc=\"Gene=";
        print "$gene";
        print "\|$cgc"              if ($cgc);
        print "; ID=$prot_id"       if ($prot_id);
        print "\";";

        print " species=\"Caenorhabditis $species\"";
        print "\n";
    }

    # Importing header lines from wormpep (see below [1] for details):
    elsif ( ($species eq 'elegans') 
             and ($input =~ /\A > (\S+)                  # $1    == $cds
                             \s+ (CE\d+)                 # $2    == $prot
                             \s+ (WBGene\d+)             # $3    == $gene
                             \s+ (locus:(\S+)            # $5    == $cgc
                             )?                          # [$4]  == [optional element]
                             \s* (\S .+? \S )?           # $6    == $blurb
                             \s* status: \w+ \s* 
                             ( (SW|TR) : (\S+)           # $9    == $uniprot
                             )?\s* ( protein_id: (\S+)   # $11   == $prot_id
                             )?                          # [$10] == [optional element]
                             /x) ) {                     # ignore whitespace in regex
        $cds     = $1;
        $prot    = $2;
        $gene    = $3;
        $cgc     = $5;
        $blurb   = $6;
        $uniprot = $9;
        $prot_id = $11;

        # CDS to mol. gene names: "T13A10.10a" -> "T13A10.10"
        if ($cds =~ /^(\S+?)[a-z]*$/) {   # use nongreedy S+ !
            $mgen = $1;
        }

        # Exporting InParanoid-ized header lines (see [3] below):
        print ">$prot";
        print " source=$source;";
        print " version=$version;";
        print " symbol=$cgc;"       if ($cgc);
        print " symbol:$mgen;"      if ((! $cgc) and ($mgen));
        print " uniprot=$uniprot;"  if ($uniprot);

        print " desc=\"Gene=";
        print "$gene\|$mgen";
        print "\|$cgc"              if ($cgc);
        print "; ID=$prot_id"       if ($prot_id);
        print "; CDS=$cds";
        print "; $blurb"            if ($blurb);
        print "\";";

        print " species=\"Caenorhabditis $species\"";
        print "\n";
    }

    elsif ( ($species eq 'remanei')
             and ($input =~ /\A > ( (.+)       # $2    == $gene
                                    .\d+)      # $1    == $prot
                                    \s+ 
                                    (\S.+\S)   # $3    == $blurb
                                    \s* 
                                    /x) ) {    # ignore whitespace in regex
        $prot  = $1;
        $gene  = $2;
        $blurb = $3;

        # Exporting InParanoid-ized header lines (see [2] below):
        print ">$prot";
        print " source=$source;";
        print " version=$version;";
        print " symbol=$gene;";

        print " desc=\"Gene=$gene";
        print "; $blurb"            if ($blurb);
        print "\";";

        print " species=\"Caenorhabditis $species\"";
        print "\n";
    }

    elsif ($input =~ /^>/) { 
        die "Failed to parse header line: $input\n";
    }
}


# [1]. Importing header lines from brigpep2gid:
# 
# Currently, this has no isoforms and thus doesn't need to pass on information
#    about isoforms in the blurb.  That might change, though, with the next
#    C. briggsae gene prediction set.
# 
# Current sample input lines:
# 
# >CBP09851       WBGene00000333|CBG14136 SW:P32810       protein_id:CAE68380.1
# >CBP18174       WBGene00000343|CBG14138 TR:Q619S2       protein_id:CAE68382.1
# >CBP10328       WBGene00000365|CBG16767
# >CBP00003       WBGene00023521|CBG00001 TR:Q629I0       protein_id:CAE57167.1

# [2]. Importing header lines from wormpep:
# 
# Input format:
#    >cosmid.number CE\d+ WBGene\d+ [locus:\S+] \
#        [maybe boilerplate] status:\S+ [SW|TR:\S+] protein_id:\S+
#    [n.b. recent proteins may be missing some of these -- e.g., lin-3 isoform below.]
#
# Sample input:
#   >C07H6.7 CE03975 WBGene00003024 locus:lin-39 \
#       C. elegans homeobox gene lin-39 status:Confirmed SW:P34684 protein_id:AAK85445.1
#   >F36H1.4d CE40413 WBGene00002992 locus:lin-3 status:Partially_confirmed

# [3]. Exporting InParanoid-ized header lines:
# 
# Output format:
#     ><MOD_identifer> source=<source_db>; version=<version>; \
#         symbol=<symbol>; uniprot=<uniprot_acc>; desc="<description>"
#
# Sample output:
#
#     >CE15373 source=WormBase; version=170; symbol=trxr-2; \ 
#         uniprot=P30635; desc="Gene=WBGene00014028|ZK637.10|trxr-2; \ 
#         ID=CAA77459.1; CDS=ZK637.10; Glutathione reductase"; \ 
#          species="Caenorhabditis elegans"


