#!/usr/bin/env perl 

# pfam_tbl2roster_nov2009.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/10/2009.
# Purpose: (LEGACY VERSION) given Pfam-A.hmm, cds2gene_Nspp.tsv, and PROTEOME.hmmscan.Pfam-A.tbl, make descending genes/PFAM roster.
# Note: uses now-obsolete table outputs from past version of HMMER 3.0.

use strict;
use warnings;

my $pfam_a     = $ARGV[0];
my $prot2gene  = $ARGV[1];
my $pfam_tbl   = $ARGV[2];

my $name       = q{};
my $acc        = q{};
my $desc       = q{};

my %name2acc   = ();
my %name2desc  = ();
my %prot2gene  = ();

my $name2gene_ref;

my %name2count = ();

my %obs_names  = ();
my %obs_genes  = ();


if ( (! $pfam_a ) or (! $prot2gene) or (! $pfam_tbl) ) { 
    die "Format: pfam_tbl2roster.pl [Pfam-A.hmm] [cds2gene_Nspp.tsv] [PROTEOME.hmmscan.Pfam-A.tbl]\n";
}

open my $PFAM_A, '<', $pfam_a or die "Can't open Pfam-A.hmm $pfam_a: $!";
while (my $input = <$PFAM_A>) { 
    chomp $input;

# Read in ordered trios of NAME to ACC and DESC; die if violated.  Sample input:
# 
# NAME  1-cysPrx_C
# ACC   PF10417.2
# DESC  C-terminal domain of 1-Cys peroxiredoxin

    if ( $input =~ /\A NAME \s+ (.+) \z /xms ) { 
        $name = $1;
        $name =~ s/\s+\z//;
        $acc  = q{};
        $desc = q{};
    }
    if ( $input =~ /\A ACC \s+ (\S+)/xms ) { 
        $acc  = $1;
    }
    if ( $input =~ /\A DESC \s+ (\S.*\S) \s* \z/xms ) { 
        $desc = $1;
        if ( (! $name) or (! $acc) or (! $desc) ) {
            die "Can't parse Pfam-A.hmm $pfam_a at: $input\n";
        }
        if ( ($name) and ($acc) and ($desc) ) { 
            $name2acc{$name} = $acc;
            $name2desc{$name} = $desc;
            $name = q{};
            $acc  = q{};
            $desc = q{};
        }
    }
}
close $PFAM_A or die "Can't close filehandle to Pfam-A.hmm $pfam_a: $!";

open my $PROT2GENE, '<', $prot2gene or die "Can't open cds2gene_Nspp.tsv file $prot2gene: $!";
while (my $input = <$PROT2GENE>) { 
    chomp $input;
    my $prot = q{};
    my $gene = q{};

# Sample input:
# 
# g10026.t1       g10026
# g10026.t2       g10026
# g10026.t3       g10026

    if ($input =~ /\A (\S+) \s+ (\S+) /xms ) { 
        $prot = $1;
        $gene = $2;
        $prot2gene{$prot} = $gene;
    }
    if ($input !~ /\A \S+ \s+ \S+ /xms ) {
        die "From file $prot2gene, can't parse input: $input\n";
    }
}
close $PROT2GENE or die "Can't close filehandle to cds2gene_Nspp.tsv file $prot2gene: $!";

open my $PFAM_TBL, '<', $pfam_tbl or die "Can't open PROTEOME.hmmscan.Pfam-A.tbl file $pfam_tbl: $!";
while (my $input = <$PFAM_TBL>) { 
    chomp $input;
    my $name = q{};
    my $prot = q{};
    my $gene = q{};

# Sample input:
# 
#   #                                         --- full sequence ----
#   # target             query                  E-value  score  bias
#   #------------------- -------------------- --------- ------ -----
#   DUF1459              g1.t1                  7.1e-08   31.9  14.4
#   Collagen             g10.t1                 3.3e-15   55.0 127.2

    if ( ( $input !~ /\A [#] /xms ) 
         and ( $input =~ /\A (\S+) \s+ (\S+) \s+ \S /xms ) ) { 

        # Capture data:
        $name = $1;
        $prot = $2;

        # Enforce reality of data:
        if (! exists $prot2gene{$prot} ) {
            die "Protein $prot has no defined gene!\n";
        }  
        if (    (! exists $name2acc{$name}  )
             or (! exists $name2desc{$name} ) ) { 
            die "PFAM-A domain $name has no defined accession number $name2acc{$name} or description $name2desc{$name}!\n";
        }

        # Store relevant data:
        $gene = $prot2gene{$prot};
        $name2gene_ref->{$name}->{$gene} = 1;
        $obs_genes{$gene} = 1;
        $obs_names{$name} = 1
    }
}
close $PFAM_TBL or die "Can't close filehandle to PROTEOME.hmmscan.Pfam-A.tbl file $pfam_tbl: $!";

foreach my $name (sort keys %{ $name2gene_ref } ) { 
    my @genes = sort keys %{ $name2gene_ref->{$name} };
    my $count = scalar @genes;
    $name2count{$name} = $count;
}

my $gene_total = scalar(keys %obs_genes);
my $pfam_total = scalar(keys %obs_names);
$gene_total    = commify($gene_total);
$pfam_total    = commify($pfam_total);

print "Total genes:          $gene_total\n";
print "Total PFAM-A domains: $pfam_total\n";
print "\n";

foreach my $name ( sort { $name2count{$b} <=> $name2count{$a} } 
                   keys %name2count ) { 
    if (    (! exists $name2count{$name} ) 
         or (! exists $name2acc{$name}   ) 
         or (! exists $name2desc{$name}  ) ) { 
        die "Can't parse PFAM domain $name in final report!\n";
    }
    my $count = $name2count{$name};
    my $acc   = $name2acc{$name};
    my $desc  = $name2desc{$name};
    print "$count\t$acc\t$name\t$desc\n";
}


sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}

