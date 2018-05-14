#!/usr/bin/env perl

# ali_uniqs2ncRNAs.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/3/2008.
# Purpose: map hits in an Ali *uniqs*.bed file to ncRNA sites in a .tsv file.

use strict;
use warnings;

unless ($#ARGV == 1) { 
    die 'Format: ./ali_uniqs2genes.pl [uniqs].bed  *.tsv ', "\n";
}

my $uniqs_file = $ARGV[0];
my $tsv_file   = $ARGV[1];
my $chromosome = 'I|II|III|IV|V|X';

# Record which nucleotides in the genome actually have 1+ uniq hits.
# Any site not recorded now, gets ignored in later GFF3 upload.  
# N.B.: This doesn't give uniqs/genelist!  Have to do that later.

my %coord2uniqs_ref = ();

open my $UNIQS, "<", "$uniqs_file" 
     or die "Can't open ", '[uniqs].bed file ', "$uniqs_file\n";

while (my $input = <$UNIQS>) { 
    chomp $input;

    # Typical input, as of 11/3/2008:
    # chrI  112081  112112  FC12873-381658  1000  +   -   -   0,0,255  1  31  0
    # chrI  112081  112112  FC12873-409950  1000  -   -   -   255,0,0  1  31  0

    if ( $input =~ /\A chr($chromosome) \s+ (\d+) \s+ (\d+) /xms ) { 
        my ($chr, $nt1, $nt2) = ($1, $2, $3);
        foreach my $i ($nt1..$nt2) { 
            $coord2uniqs_ref{$chr}->{$i} = 1;
        }
    }
}

close $UNIQS;

# For each nt with a uniq, list any gene(s) with an exon covering it.

my %coord2genes_ref = ();

open my $TSV, "<", "$tsv_file" 
    or die "Can't open TSV file $tsv_file\n";

while (my $input = <$TSV>) { 
    chomp $input;

    # Typical (tab-delimited) input:
    # WBGene00003039  lin-58  F56A12.3        miRNA   V       14364450        14364472

    if ( $input =~ /  \A 
                          (WBGene\d+) 
                      \s+ (\S+) 
                      \s+ (\S+) 
                      \s+ \S+ 
                      \s+ (\S+)  # ($chromosome) 
                      \s+ (\d+) 
                      \s+ (\d+)
                   /xms ) {
        my ($wbgene, $cgc, $cds, $chr, $nt1, $nt2);
        ($wbgene, $cgc, $cds, $chr, $nt1, $nt2) = ($1, $2, $3, $4, $5, $6);
        my $gene = $wbgene . '|' . $cds . '|' . $cgc;
        foreach my $i ($nt1..$nt2) { 
            if (exists $coord2uniqs_ref{$chr}->{$i} ) { 
                $coord2genes_ref{$chr}->{$i}->{$gene} = 1;  
            }
            # I.e.: for nt $i of chromosome $chr, 
            #   if it got recorded as having a signal!
            #   and, if $gene has an exon, store as hash key.
            #   N.B.: >1 genes can be stored this way.
        }
    }
}

close $TSV;

# Now, reread uniqs, but this time, assign a count of each uniq to 
#    any genelist it touches (ideally just one list with one gene),
#    but count overlaps as belonging to both classes.

my %gl2score     = ();

open $UNIQS, "<", "$uniqs_file"
     or die "Can't open *uniqs*.bed file $uniqs_file\n";

while (my $input = <$UNIQS>) {
    chomp $input;
    
    # Typical input, as of 11/3/2008:
    # chrI  112081  112112  FC12873-381658  1000  +   -   -   0,0,255  1  31  0
    # chrI  112081  112112  FC12873-409950  1000  -   -   -   255,0,0  1  31  0
 
    if ( $input =~ /\A chr($chromosome) \s+ (\d+) \s+ (\d+) /xms ) {
        my $chr;
        my $nt1;
        my $nt2;
        ($chr, $nt1, $nt2) = ($1, $2, $3);
        my %glists_seen = ();
        foreach my $i ($nt1..$nt2) { 
            my @genes = sort keys %{ $coord2genes_ref{$chr}->{$i} };
            my $genelist = join '; ', @genes;
            $glists_seen{$genelist} = 1;
        }
        foreach my $glist (sort keys %glists_seen) { 
            $gl2score{$glist}++;
        }
    }
}
 
close $UNIQS;

# List the (genelist) keys of %gl2score by descending signal strength.
# N.B.: for some reason, ghost counts are getting in with a ' ' key;
#    'grep { $_ =~ /\S/ }' weeds these out.

my @readouts = sort { $gl2score{$b} <=>  $gl2score{$a} } 
                   grep { $_ =~ /\S/ } keys %gl2score;

# Print the summary.

foreach my $glist (@readouts) { 
    print "$glist\t$gl2score{$glist}\n";
}

