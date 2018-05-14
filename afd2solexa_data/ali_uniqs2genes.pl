#!/usr/bin/perl

# ali_uniqs2genes.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/18/2008.
# Purpose: map hits in an Ali *uniqs*.bed file to genes from a GFF file.

use strict;
use warnings;

unless ($#ARGV == 1) { 
    die 'Format: ./ali_uniqs2genes.pl *uniqs*.bed  *.gff3 ', "\n";
}

my $uniqs_file      = $ARGV[0];
my $gff3_file       = $ARGV[1];
my $chromosome      = 'I|II|III|IV|V|X';

# Record which nucleotides in the genome actually have 1+ uniq hits.
# Any site not recorded now, gets ignored in later GFF3 upload.  
# N.B.: This doesn't give uniqs/genelist!  Have to do that later.

my %coord2uniqs_ref = ();

open my $UNIQS, "<", "$uniqs_file" 
     or die "Can't open ", '*uniqs*.bed file ', "$uniqs_file\n";

while (my $input = <$UNIQS>) { 
    chomp $input;

    # Typical input:
    # chrIII 4753131 4753155 read 0 + - - 0,0,255

    if ( $input =~ /\A chr($chromosome) \s+ (\d+) \s+ (\d+) /xms ) { 
        my $chr;
        my $nt1;
        my $nt2;
        ($chr, $nt1, $nt2) = ($1, $2, $3);
        foreach my $i ($nt1..$nt2) { 
            $coord2uniqs_ref{$chr}->{$i} = 1;
        }
    }
}

close $UNIQS;

# For each nt with a uniq, list any gene(s) with an exon covering it.

my %coord2genes_ref = ();

open my $GFF3, "<", "$gff3_file" 
    or die "Can't open GFF3 file $gff3_file\n";

while (my $input = <$GFF3>) { 
    chomp $input;

    # Typical (tab-delimited) input:
    # IV  Coding_transcript  exon  11054282  11054340  .  + . \\ 
    #     Parent=Gene:WBGene00002992

    if ( $input =~ /  \A 
                      ($chromosome)     \t 
                      Coding_transcript \t 
                      exon              \t 
                      (\d+)             \t
                      (\d+)             \t
                      (?: [^\t]* \t ){3}
                      Parent=Gene: 
                      (WBGene\d+) 
                   /xms ) {
        my ($chr, $nt1, $nt2, $gene);
        ($chr, $nt1, $nt2, $gene) = ($1, $2, $3, $4);
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

close $GFF3;

# Now, reread uniqs, but this time, assign a count of each uniq to 
#    any genelist it touches (ideally just one list with one gene),
#    but count overlaps as belonging to both classes.

my %gl2score     = ();

open $UNIQS, "<", "$uniqs_file"
     or die "Can't open *uniqs*.bed file $uniqs_file\n";

while (my $input = <$UNIQS>) {
    chomp $input;
    
    # Typical input:
    # chrIII 4753131 4753155 read 0 + - - 0,0,255
 
    if ( $input =~ /\A chr($chromosome) \s+ (\d+) \s+ (\d+) /xms ) {
        my $chr;
        my $nt1;
        my $nt2;
        ($chr, $nt1, $nt2) = ($1, $2, $3);
        my %glists_seen = ();
        foreach my $i ($nt1..$nt2) { 
            my @genes = sort keys %{ $coord2genes_ref{$chr}->{$i} };
            my $genelist = join '|', @genes;
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

