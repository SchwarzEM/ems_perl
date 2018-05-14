#!/usr/bin/env perl

# minim_gff2_to_gff3.pl -- Erich Schwarz <emsch@its.caltech.edu>, 5/25/2009.
# Purpose: map minimal GFF2 to minimal GFF3.

use strict;
use warnings;
use Getopt::Long;

# Set this to '1' to see what genes are being skipped:
my $DEBUG;
my $HELP;
# Note that this script, as written, skips ncRNAs, psuedogenes, and uncloned loci.

# The default is what I *think* is correct GFF3 syntax:
my $exon_name = 'coding_region_of_exon';
# But to pass WormBase 'validation', I may need to set this to 'CDS' instead.

GetOptions ( "help"        => \$HELP,
             "debug"       => \$DEBUG,
             "exon_name=s" => \$exon_name, );

if ($HELP) { 
    die "Format: minim_gff2_to_gff3.pl [--help] [--debug] [--exon_name TEXT] input.gff2 [to STDOUT]\n";
}

# Die loudly if given a bad exon name:
if ( $exon_name !~ /\A \w+ \z/xms ) { 
    die "Will not accept \"$exon_name\" as a term",
        " for a \"CDS\", i.e., a coding region of an exon!\n",
        ;
}

# Link CDSes to transcripts: note that one CDS can have >1 Transcripts...
my $cds2tx_ref;   # Don't define as 'q{}' if it's a hashref.

# Link transcripts to genes.
my %tx2gene     = ();

# Be able to say if a gene is transcribed (and, if need be, 
#     if it has a particular tx, although that capacity doesn't 
#     get used in this script).
my $gene2tx_ref;

my $input_gff = $ARGV[0];

# Do first round just to link CDSes to transcripts.
# 
# Note that CDSes are themselves pretty useless in GFF3 since they're not meant
#    in anything like the sloppy way that they're used in WormBase GFF2s.
# But, to keep an unbroken chain of identities, we have to extract the 
#    CDS->Transcript links from the GFF2.

open my $GFF2, '<', $input_gff
    or die "Can't open input GFF2 $input_gff: $!";
while (my $input = <$GFF2>) {

# Sample input:
# Cbre_Contig86   Coding_transcript       coding_exon     217809  217958  .       -       0       Transcript "CBN06464" ; CDS "CBN06464"

    chomp $input;
    my $cds1 = q{};
    my $tx1  = q{};
    if ( $input =~ /\A \S+ \t 
                     Coding_transcript \t 
                     coding_exon \t /xms ) { 
        if ( $input =~ / Transcript \s+ \" ( [^\"]+? ) \" /xms ) {
            $tx1 = $1;
        }
        if ( $input =~ / CDS \s \" ( [^\"]+? ) \" /xms ) {
            $cds1 = $1;
        }
        if ( $tx1 and $cds1 ) {
            # Unfortunately, WormBase can have 2+ Transcripts encoding 1 CDS:
            $cds2tx_ref->{$cds1}->{$tx1} = 1;
        }
    }
}
close $GFF2 or die "Can't close input GFF2 $input_gff: $!";

# Do second round to get CDS->Transcript names linked to genes!

open $GFF2, '<', $input_gff 
    or die "Can't open input GFF2 $input_gff: $!";
while (my $input = <$GFF2>) { 

# Sample input:
# Cbre_Contig86   curated CDS     215137  217958  .       -       .       CDS "CBN06464" ; WormPep "CN:CN27269" ;  Status "Predicted" ;  Gene "WBGene00145189" ;

    chomp $input;
    my $gene = q{};
    my $cds2 = q{};
    my $tx2  = q{};
    if ( $input =~ /\A \S+ \t curated \t CDS \t /xms ) { 
        if ( $input =~ / Gene \s \" ( WBGene\d+ ) \" /xms ) { 
            $gene = $1;
        } 
        if ( $input =~ / CDS \s \" ( \S+ ) \" /xms ) {
            $cds2 = $1;
        }
        if ( ( $gene and $cds2 ) and ( exists $cds2tx_ref->{$cds2} ) ) { 
            foreach my $listed_tx ( sort keys %{ $cds2tx_ref->{$cds2} } ) { 
                $gene2tx_ref->{$gene}->{$listed_tx} = 1; 
                $tx2gene{$listed_tx} = $gene;
            }
        }
    }
}

close $GFF2 or die "Can't close input GFF2 $input_gff: $!";

# At last, do third round to convert+print easy lines, generate+print hard lines.
open $GFF2, '<', $input_gff
    or die "Can't open input GFF2 $input_gff: $!";

# Begin with header which I hope is valid...
print "##gff-version 3\n";
         
while (my $input = <$GFF2>) {  
    chomp $input;

    # Reformat 'commented' sequence data to real GFF3 sequence lines.
    # Sample input: ##sequence-region Cbre_Contig0 1 4147112

    if ( $input =~ / \A [#][#]sequence-region \s+ \S+ \s \d+ \s+ \d+ /xms ) { 
        print "$input\n";
    }

    # For now, I'm skipping all other comment lines -- simplest way to avoid errors.

    # Map 'coding_exon' to $exon_name -- default value, 'coding_region_of_exon':
    # 
    # Sample input:
    # Cbre_Contig86   Coding_transcript       coding_exon     217809  217958  .       -       0       Transcript "CBN06464" ; CDS "CBN06464"

    if ( $input =~ / \A 
                     ( \S+ \t 
                       Coding_transcript \t 
                       coding_exon \t 
                       \d+ \t \d+ \t \. \t \S \t \d \t )  # $1
                     Transcript \s+ \" ([^\"]+) \" [^\t]* # $2
                     \z /xms ) { 
        my $line1 = $1;
        my $tx3   = $2;

        # It is this, and not 'CDS', which is the correct mapping into GFF3!
        $line1 =~ s/coding_exon/$exon_name/;

        # Note: $line1 has a terminal \t already; don't add another!
        print "$line1",

              # Note: I am not even trying to give the individual exons a name.
              # I am also totally avoiding the questionable practice of considering 
              #    a string of non-contiguous sequence blocks to be aggregated into 
              #    a 'CDS' just because it works groovy that way in WormBase GFF2.
              # Instead, I am just giving the exons a single parent.

              "Parent=$tx3\n",
              ;
    }

    # Generate mRNA lines:
    # 
    # Sample input:
    # Cjap_Contig63   Coding_transcript       protein_coding_primary_transcript       143286  153269  .       +       .       Transcript "CJA04690"

    # This only works if %tx2gene has already been populated.

    if ( $input =~ / \A ( \S+ \t Coding_transcript \t protein_coding_primary_transcript \t .+ \t )
                        Transcript \s+ \" ( [^\"]+? ) \" [^\t]* \z /xms ) { 
        my $line2 = $1;
        my $tx4  = $2;
        $line2 =~ s/protein_coding_primary_transcript/mRNA/;
        # Again, $line2 already ends in '\t':
        print "$line2", "ID=$tx4";  # Not "ID=Transcript:$tx4", ";Parent=Gene:...".
        print ";Parent=$tx2gene{$tx4}" if ( exists $tx2gene{$tx4} );
        if (! exists $tx2gene{$tx4} ) { 
            warn "Transcript $tx4 has no known Gene parent!\n";
        }
        print "\n";
    }

    # Sample input:
    # chrII   gene    gene    716619  731726  .       -       .       Gene "WBGene00024225" ; Locus "Cbr-let-23" ; Locus "Cbr-let-23"
    # 
    # Only do this mapping *if* an observed transcript exists for the gene!

    if ( $input =~ / \A ( \S+ \t gene \t gene \t .+ \t )  # $1
                     Gene \s+ \" ([^\"]+) \" [^\t]* \z    # $2
                   /xms ) {
        my $line3 = $1;
        my $gene1 = $2;
        $line3 =~ s/gene\tgene/Coding_transcript\tgene/;
        if (exists $gene2tx_ref->{$gene1} ) { 
            # And $line3 already ends w/ '\t':
            print "$line3", "ID=$gene1\n";     # Not "ID=Gene:".
        }
        if (! exists $gene2tx_ref->{$gene1} and $DEBUG ) {
            warn "Tried to print line for Gene $gene1 that, apparently, has no transcript!\n";
        }
    }
}

close $GFF2 or die "Can't close input GFF2 $input_gff: $!";

