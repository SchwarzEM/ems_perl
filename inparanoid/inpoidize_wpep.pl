#!/usr/bin/perl

# inpoidize_wpep.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/5/2007
# Purpose: convert standard wormpep headers into InParanoid-preferred headers.

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

unless ($#ARGV == 2) { 
    die "Format: ./inpoidize_wpep.pl [fasta] [source DB] [version]\n";
}
my $version = pop @ARGV;
my $source  = pop @ARGV;

while ($input = <>) { 
    # Unless it's a header line, just pass on unchanged:
    if ($input !~ /^>/) { 
        print $input;
    }

    # Importing header lines from wormpep (see below [1] for details):
    elsif ($input =~ /\A > (\S+)                             # $1    == $cds
                      \s+ (CE\d+)                            # $2    == $prot
                      \s+ (WBGene\d+)                        # $3    == $gene
                      \s+ (locus:(\S+)                       # $5    == $cgc
                      )?                                     # [$4]  == [optional element]
                      \s* (\S .+? \S )?                      # $6    == $blurb
                      \s* status: \w+ \s* ( (SW|TR) : (\S+)  # $9    == $uniprot
                      )?\s* ( protein_id: (\S+)              # $11   == $prot_id
                      )?                                     # [$10] == [optional element]
                      /x) {                                  # ignore whitespace
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

        # Exporting InParanoid-ized header lines (see [2] below):
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
        print "\"\n";
    }
    elsif ($input =~ /^>/) { 
        die "Failed to parse header line: $input\n";
    }
}


# [1]. Importing header lines from wormpep:
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

# [2]. Exporting InParanoid-ized header lines:
# 
# Output format:
#     ><MOD_identifer> source=<source_db>; version=<version>; \
#         symbol=<symbol>; uniprot=<uniprot_acc>; desc="<description>"
#
# Sample output:
#
#     >CE15373 source=WormBase; version=170; symbol=trxr-2; \
#         uniprot=P30635; desc="Gene=WBGene00014028|ZK637.10|trxr-2; \
#         ID=CAA77459.1; CDS=ZK637.10; Glutathione reductase"

