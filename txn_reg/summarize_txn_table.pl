#!/usr/bin/perl

# summarize_txn_table.pl -- Erich Schwarz <emsch@its.caltech.edu>, 9/4/2007.
# Purpose: convert txn_fac1_v17.txt to a list of nonredundant refs. and putative TFs.

# One example of usage:
#   ./summarize_txn_table.pl   \
#       /home/schwarz/taygeta/science2/wormbase/Ongoing/wpa_xref.txt                 \
#       /home/schwarz/taygeta/science/ali_work/macro_sets_2007/ws170_gene_names.txt  \
#       /home/schwarz/taygeta/science/txn_factor_data/txn_fac1_v17.txt > test.out

use strict;
use warnings;

my %pmid2wbpaper  = ();
my %seen_pmid     = ();
my %tfname2geneid = ();
my %seen_tf       = ();

while (my $input = <>) { 
    chomp $input;
    if ($input =~ / \[ (WBPaper\d+) .+? (pmid\d+) \] /xms) { 
        $pmid2wbpaper{$2} = $1;
    }
    elsif ($input =~ /\A WBGene\d+ \s+ [^\.\s]+ \. /xms) { 
        my $geneid = $input;
        $geneid =~ s/\t/|/g;
        if ($input =~ /\A
                       WBGene\d+
                       \s+
                       ([^\.\s]+ \. [^\.\s]+)
                       \s+
                       ([^-\s]+ - [^-\s]+) \z 
                       /xms) {
            $tfname2geneid{$1} = $geneid;
            $tfname2geneid{$2} = $geneid;
        }
        elsif ($input =~ /\A 
                       WBGene\d+ 
                       \s+ 
                       (^[\.\s]+ \. ^[\.\s]+)
                       \s* \z 
                       /xms) { 
            $tfname2geneid{$1} = $geneid;
        }
    }
    else { 
        if ($input =~ /\A (\S+(\.|-)\S+) /xms) { 
            $seen_tf{$1} = 1;
        }
        if ($input =~ /PMID:\s+(\d+)/) { 
            $seen_pmid{"pmid$1"} = 1;
        }
        if ($input =~ /(pmid\d+)/) {
            $seen_pmid{"$1"} = 1;
        }
    } 
}

foreach my $tf (sort keys %seen_tf) { 
    if ($tfname2geneid{$tf}) { 
        print "$tfname2geneid{$tf}\n";
    }
    else { 
        print "UNRECOGNIZED: $tf\n";
    }
}

foreach my $pmid (sort keys %seen_pmid) { 
    if ($pmid2wbpaper{$pmid}) { 
        print "$pmid\t$pmid2wbpaper{$pmid}\n";
    }
    else { 
        print "$pmid\n"
    }
}

