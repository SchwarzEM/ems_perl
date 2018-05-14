#!/usr/bin/env perl

# genes2eles.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/3/2011.
# Purpose: given C. elegans gene data from a GFF, and element coordinates, annotate the elements.

use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils qw(uniq);
use Roman;   # To get 'arabic' function, needed for sorting names in sub namevalue.

my %opts = ();

my @chromosomes = ();
my %OK_chr      = ();

GetOptions ( 'chromosome=s{,}' => \@chromosomes,
             'elements=s'      => \$opts{'element_table'},
             'fasta=s'         => \$opts{'element_fasta'},
             'best_mot=s'      => \$opts{'best_mot'},
             'names=s'         => \$opts{'gene_names'},
             'protein=s'       => \$opts{'protein_coding'},
             'modencode=s'     => \$opts{'modENCODE'},
             'genome=s'        => \$opts{'genome_GFF'},
             'help'            => \$opts{'help'},         );

if (    $opts{'help'} 
     or (! $opts{'element_table'}  ) 
     or (! $opts{'element_fasta'}  )
     or (! $opts{'gene_names'}     )
     or (! $opts{'protein_coding'} )
     or (! $opts{'modENCODE'}      )  
     or (! $opts{'genome_GFF'}     )
   ) { 
    die "\n",
        "Format: genes2eles.pl\n",
        "          --chromosomes|-c  [List of 1 or more chromosomes to annotate; default is to annotate all]\n",
        "          --elements|-e     [Table of coordinates of elements]\n",
        "          --fasta|-f        [FASTA of element sequences]\n",
        "          --best_mot|-b     [Table of which motifs have best-scoring motif hits]\n",
        "          --names|-n        [TSV of WBGene, sequence, and CGC gene names]\n",
        "          --protein|-p      [file with WBGene IDs of protein-coding genes -- can be wormpep, or just list]\n",
        "          --modencode|-m    [list of 7,237 WBGenes identified as new ncRNAs by modENCODE by WS220]\n",
        "          --genome|-g     [GFF with coords. of gene boundaries, exons/introns, and TSSes]\n",
        "          --help          [Print this message and quit]\n",
        "\n",
        ;
}

# These get populated in first reading of GFF, which also enforces upload of lengths.
my %chrs = ();

my $element     = q{};
my $gene        = q{};
my $tx          = q{};
my $chr         = q{};
my $chr_len     = q{};
my $ex_int      = q{};   # Exon or intron.
my $nt1         = q{};
my $nt2         = q{};
my $ori         = q{};
my $wb_geneid   = q{};
my $hum_geneid  = q{};
my $tx_start_nt = q{};
my $bad_gene    = q{};

my %gene_status = ();
my %pseudo      = ();

my $info_ref;

# Map each transcript to a unique WBGene parent, and each WBGene name to a common name; simple hash works.
my %tx2gene   = ();
my %gene2name = ();

# Map each gene to one, or more, transcripts; thus, need hashref.
my $gene2tx_ref;

# Optionally, annotate only one or a few chromosomes (useful for limiting huge RAM requirement):
if (@chromosomes) { 
    foreach my $listed_chr (@chromosomes) { 
        $OK_chr{$listed_chr} = 1;
    }
}

# Before collecting and processing annotation data, first record data for elements, 
#     so that (subsequently) only those nucleotides in the genome belonging to an element actually get recorded, 
#     thus greatly reducing load on RAM!

open my $TABLE, '<', $opts{'element_table'} or die "Can't open element-table file $opts{'element_table'}: $!";
while (my $input = <$TABLE>) {
    chomp $input;

# Example of input lines:
# I       ce6_ct_v6step13 exon    3987    4011    0.000000        +       .       gene_id "chrI.1"; transcript_id "chrI.1";
# III   pmid18981268    MUSSA   7530656 7531620 .       .       .       Cons_block "N1.1"

    if ( $input =~ / \A 
                     (\S+)                # chrI-X
                     \t [^\t]* \t [^\t]* 
                     \t (\d+)             # start nt
                     \t (\d+)             # end nt
                     \t [^\t]* 
                     \t ([+]|[.]|[-])     # orientation -- should be '+' or '.'
                     \t [^\t]* 
                     \t (?: gene_id | Cons_block ) \s+ \" ([^\"\s]+) \"  # element name
                      /xms ) { 
        $chr     = $1;
        $nt1     = $2;
        $nt2     = $3;
        $ori     = $4;
        $element = $5;

        # Only go on to the following actions if $chr is free to be anything, or does belong to a specified OK list:
        if ( (! @chromosomes) or ( exists $OK_chr{$chr} ) ) {
            # Enforce reliability of ascending coordinates and '+' orientation:
            if ( ( $nt1 > $nt2 ) or ( $ori eq '-' ) ) { 
                die "Elements should uniformly be in '+' or '.', not '-', orientation, and only with ascending coordinates: $input\n";
            }
            # Convert '.' to '+'.
            if ($ori eq q{.}) { 
                $ori = '+';
            }
            $info_ref->{'element'}->{$element}->{'chr'}   = $chr;
            $info_ref->{'element'}->{$element}->{'start'} = $nt1;
            $info_ref->{'element'}->{$element}->{'end'}   = $nt2;
            foreach my $nt ($nt1..$nt2) { 
                $info_ref->{'in_an_element'}->{'chr'}->{$chr}->{'nt'}->{$nt} = 1;
            }
        }
    }
    else { 
        if ( $input !~ / \A \# /xms ) { 
            die "Can't parse input line from element-table file $opts{'element_table'}: $input\n";
        }
    }
}
close $TABLE or die "Can't close filehandle to element-table file $opts{'element_table'}: $!";

# Get sequence characteristics for each element:
open my $ELE_FASTA, '<', $opts{'element_fasta'} or die "Can't open FASTA file of sequences of elements: $!";

# Set to zero-bytes:
$element = q{};

# Read in the data:
while (my $input = <$ELE_FASTA>) { 
    chomp $input;
    if ( $input =~ /> (\S+) /xms ) { 
        $element = $1;
        if (! exists $info_ref->{'element'}->{$element} ) { 
            # In essence, refuse to 'see' elements that haven't been entered above (maybe because of a chromosomal restriction).
            $element = q{};
        }
    }
    elsif ($element and ( $input =~ /\A ([ACGTacgt]+) \z /xms ) ) { 
        my $ele_seq = $1;
        $info_ref->{'element'}->{$element}->{'sequence'} .= $ele_seq;
    }
    else {
        # It's safe to ignore nonsense lines, as long as there's no element they're annotating! 
        if ($element) { 
            die "Can't parse input line from $opts{'element_fasta'}: $input\n";
        }
    }
}
close $ELE_FASTA or die "Can't close FASTA file of sequences of elements: $!";

# Find which elements, if any, have best-scoring motif hits:
open my $BEST_MOT, '<', $opts{'best_mot'} or die "Can't open table of elements with best-scoring motif hits $opts{'best_mot'}!\n";
while (my $input = <$BEST_MOT>) { 
    chomp $input;
    my $best_hit_motif    = q{};
    my $best_hit_text     = q{};
    my @best_hit_elements = ();

# Sample inputs:
# 2-24	3 elements: [q = 6.98e-06, 3 elements: chrI.3594, chrII.9001, chrIII.6083]
# 2-26	1 element: [q = 4.05e-06, 1 element: chrIV.16093]

    if ( $input =~ / \A
                     (\S+)
                     \s+
                     .+
                     (?: element | elements) : 
                     (.+?)
                     \] /xms ) {
        $best_hit_motif    = $1;
        $best_hit_text     = $2;
        $best_hit_text     =~ s/,//g;
        @best_hit_elements = split /\s+/, $best_hit_text;
        foreach my $best_hit_element (@best_hit_elements) { 
            if ( exists $info_ref->{'element'}->{$best_hit_element} ) { 
                $info_ref->{'element'}->{$best_hit_element}->{'best_motifs'}->{$best_hit_motif} = 1;
            }
        }
    }
    else { 
        die "Can't parse input from $opts{'best_mot'}: $input!\n";
    }
}
close $BEST_MOT or die "Can't close filehandle of of elements with best-scoring motif hits $opts{'best_mot'}!\n";

foreach my $element1 ( sort keys %{ $info_ref->{'element'} } ) { 
    if ( exists $info_ref->{'element'}->{$element1}->{'best_motifs'} ) { 
        my @best_motifs = sort keys %{ $info_ref->{'element'}->{$element1}->{'best_motifs'} };
        $info_ref->{'element'}->{$element1}->{'best_motif_desc'} 
            = join '; ', @best_motifs;
    }
    else { 
        $info_ref->{'element'}->{$element1}->{'best_motif_desc'} = q{};
    }
}

# Enforce completeness of element data:
foreach my $element1 ( sort keys %{ $info_ref->{'element'} } ) { 
    if (    (! exists $info_ref->{'element'}->{$element1}->{'chr'}      ) 
         or (! exists $info_ref->{'element'}->{$element1}->{'sequence'} ) ) { 
        die "Element $element1 does not have both an assigned chromosome and a defined sequence!\n";
    }
}

# Work up a human-readable description of each element's sequence:
foreach my $element2 ( sort keys %{ $info_ref->{'element'} } ) { 
    my $ele_seq = $info_ref->{'element'}->{$element2}->{'sequence'};
    my $ele_seq_length = length($ele_seq);
    # Partly this is a sanity check, partly an enforcement of no division by zero:
    if ( $ele_seq_length < 1 ) { 
        die "Insane $element2 sequence is less than 1 residue long!\n";
    }
    # Use truly turgid Perl trick to count the number of c, g, C, or G in the sequence line:
    my $ele_cg_count    = ( $ele_seq =~ tr/[cgCG]/[cgCG]/ );
    my $ele_cg_fraction = ( $ele_cg_count / $ele_seq_length );
    $ele_cg_fraction    = sprintf("%.2f", $ele_cg_fraction);
    $info_ref->{'element'}->{$element2}->{'seq_length'} = "$ele_seq_length nt";
    $info_ref->{'element'}->{$element2}->{'seq_CG'} 
        = $ele_cg_fraction
          . q{ (} 
          . $ele_cg_count 
          . q{/} 
          . $ele_seq_length 
          . q{) CG}
          ;
}

open my $NAMES, '<', "$opts{'gene_names'}" or die "Can't open TSV of gene names $opts{'gene_names'}: $!";
while (my $input = <$NAMES>) { 
    chomp $input;
    $wb_geneid  = q{};
    $hum_geneid = q{};
    if ( $input !~ /\A WBGene\d+ \t [^\t]+ \t [^\t]* \z /xms ) { 
        die "Unacceptable format in line of gene-name TSV: $input\n";
    }
    elsif ( $input =~ /\A (WBGene\d+) \t [^\t]+ \t ([^\t]+) \z /xms ) {
        $wb_geneid  = $1;
        $hum_geneid = $2;
        if ( $hum_geneid !~ /\-/xms ) { 
            die "$hum_geneid is not a CGC name for $wb_geneid in: $input\n";
        }
    }
    elsif ( $input =~ /\A (WBGene\d+) \t ([^\t]+) \t \z /xms ) {
        $wb_geneid  = $1;
        $hum_geneid = $2;
        if ( ( $hum_geneid =~ /\-/xms ) or ( $hum_geneid !~ /\./xms ) ) {
            die "$hum_geneid is not a sequence name for $wb_geneid in: $input\n";
        }
    }
    else { 
        die "Can't parse line of gene-name TSV: $input\n";
    }
    if ( ( $wb_geneid !~ /\S/xms ) or ( $hum_geneid !~ /\S/xms ) ) {
        die "Both WBGene ID '$wb_geneid' and common name '$hum_geneid' need at least one non-space character!\n";
    }
    # Put in both types of ID to make grepping lines easier (WBGene IDs are reliable for grep; human-readable aren't!).
    $gene2name{$wb_geneid} = $wb_geneid . q{|} . $hum_geneid;

    # Default assignment of gene status is to non-modENCODE, non-protein-coding, non-pseudogene -- i.e., to ordinary ncRNA:
    $gene_status{$wb_geneid} = 'ncRNA';
}
close $NAMES or die "Can't close filehandle to TSV of gene names $opts{'gene_names'}: $!";

open my $PROTCODING, '<', "$opts{'protein_coding'}" or die "Can't open file with WBGene IDs of protein-coding genes: $!";  
while (my $input = <$PROTCODING>) {
    chomp $input;
    if ( $input =~ / (WBGene\d+) /xms ) {
        $gene = $1;
        if (! exists $gene2name{$gene} ) {
            die "Unrecognized putative protein-coding gene $gene in file $opts{'protein_coding'}, in text: $input\n";
        }
        $gene_status{$gene} = 'prot_coding';
    }
}
close $PROTCODING or die "Can't close filehandle to file with WBGene IDs of protein-coding genes: $!";

open my $MODENCODE, '<', "$opts{'modENCODE'}" or die "Can't open list of WBGenes for 7K new modENCODE ncRNAs: $!";
while (my $input = <$MODENCODE>) { 
    chomp $input;
    if ( $input =~ /\A WBGene\d+ \z/xms ) { 
        if (! exists $gene2name{$input} ) { 
            die "Unrecognized putative modENCODE 7K ncRNA gene $input in file $opts{'modENCODE'}, in text: $input\n";
        }
        $gene_status{$input} = 'modENCODE';
    }
    else { 
        die "Can't parse input line from $opts{'modENCODE'}: $input\n";
    }
}
close $MODENCODE or die "Can't close filehandle to list of WBGenes for 7K new modENCODE ncRNAs: $!";

# Multiple rounds of reading the genome GFF are needed.
# In the first round of reading:
#     find which chromosomes/contigs we actually can name;
#     find extents of them in nucleotides.

open my $GENOME, '<', "$opts{'genome_GFF'}" or die "Can't open genome GFF file $opts{'genome_GFF'}: $!";

CHR_ID_AND_LENGTH: 
while (my $input = <$GENOME>) {
    chomp $input;

    # Initialize to null values, for each new text line to be evaluated:
    $chr     = q{};
    $chr_len = q{};

    # Example of input lines for getting lengths of contigs/chromosomes:
    #
    #  ##sequence-region CHROMOSOME_MtDNA 1 13794
    #  CHROMOSOME_MtDNA    .       Sequence        1       13794   .       +       .       Sequence "MTCE"

    if ( $input =~ / \A 
                     [#]{2}
                     sequence-region 
                     \s+
                     (\S+) 
                     \s+ 
                     1 
                     \s+ 
                     (\d+) 
                     \s* \z /xms ) {
        $chr     = $1;
        $chr_len = $2;
     }

     # Input:  I	.	Sequence	1	15072423	.	+	.	Sequence "I"
     elsif ( $input =~ / \A
                      (\S+)
                      \t
                      \. 
                      \t 
                      Sequence 
                      \t 
                      1 
                      \t 
                      (\d+)
                      \t \. \t [+] \t \. \t Sequence \s+ \" [^\"]+ \" \s* \z /xms ) { 

        $chr     = $1;
        $chr_len = $2;
    }
    else { 
        next CHR_ID_AND_LENGTH; # Go on to next line of input if nothing useful found.
    }

    # Only go on to the following actions if $chr *was* set to a nonzero value and is OK:
    if ( (@chromosomes) and (! exists $OK_chr{$chr} ) ) { 
        next CHR_ID_AND_LENGTH;
    }

    # Block bad pattern-matching.
    if ( $chr =~ /\A \s* \z/xms ) {
        die "Chromosome $chr has no recognizable characters in input: $input\n";
    }   

    # We don't want two or more such lines in a file.
    if (exists $chrs{$chr}) { 
        die "On chromosome $chr, redundant specification of chromosome/contig length: $input\n";
    }

    # One-time entry of chr's existence and length, if $chr is actual text and is OK.
    if ( $chr =~ /\S/xms ) { 
        $chrs{$chr} = 1;
        $info_ref->{'genome'}->{$chr}->{'length'} = $chr_len;
    }
}

# In the second round of reading, identify Pseudogenes for later censorship:

open $GENOME, '<', "$opts{'genome_GFF'}" or die "Can't open genome GFF file $opts{'genome_GFF'}: $!";
while (my $input = <$GENOME>) {
    chomp $input;
    # III	Pseudogene	Pseudogene	1743240	1743311	.	-	.	Pseudogene "Y22D7AR.t1"
    if ( $input =~ / Pseudogene [ ]+ \" ([^\"]+) \" /xms ) { 
        $bad_gene = $1;
        $gene_status{$bad_gene} = 'pseudogene';
    }
}
close $GENOME or die "Can't close filehandle to genome GFF file $opts{'genome_GFF'}: $!";

# In the third round of reading, map genes (with WBGene\d+ IDs) to transcripts, and get boundaries of genes.
open $GENOME, '<', "$opts{'genome_GFF'}" or die "Can't open genome GFF file $opts{'genome_GFF'}: $!";

LOOP1:
while (my $input = <$GENOME>) {
    chomp $input;

    # Initialize for each text line to be evaluated:
    $chr = q{};
    $nt1 = q{};
    $nt2 = q{};
    $ori = q{};
    $tx  = q{};
    $gene = q{};

# Example of input line for linking transcripts ('CDS') to genes:
#
# II	curated	CDS	9199129	9206878	.	+	.	
# 	CDS "ZK1067.1b" ; WormPep "CE:CE42891" ;  Locus "let-23" ;  Status "Confirmed" ;  Gene "WBGene00002299" ; 
# 
# I	curated_miRNA	miRNA_primary_transcript	1738637	1738735	.	+	.	Transcript "Y71G12B.29" ; Gene "WBGene00003278" ; 

    if ( $input =~ / \A 
                     (\S+) 
                     \t .+ \t
                     # This entire block of text should be "\t"-free, hence the use of "[^\t]*" and "[ ]+".
                     (?:CDS|Transcript) [ ]+ \" ([^\"\s]+) \" 
                     [^\t]* [ ]+ ; [ ]+ 
                     Gene [ ]+ \" ([^\"\s]+) \" \s* ; 
                   /xms ) { 
        $chr  = $1;
        $tx   = $2;
        $gene = $3;

        # If this $tx belongs to a pseudogene, ban it from all further consideration:
        if ( $gene_status{$gene} eq 'pseudogene' ) { 
            $pseudo{$tx} = 1;
            next LOOP1;
        }

        if ( (@chromosomes) and (! exists $OK_chr{$chr} ) ) {
            next LOOP1;
        }

        # WormBase now has some transcripts with names like "Y74C9A.2.3", and "F53G12.5a.1", which are gumming up the works.
        # Since I have no real reason to care about specific identity of the transcripts
        #     -- I end up mapping their starts to gene names anyway! --
        #     there is no good reason to keep the second suffix of these sorts of terms.
        # So, get rid of them before doing anything else.

        if ( $tx =~ /\A ([^\.\s]+ \. \d+ [a-z]*) \. \d+ /xms ) { 
            $tx = $1;
        }

        # Each tx-to-gene mapping is unique:
        $tx2gene{$tx} = $gene;

        # But a single gene can have multiple transcripts:
        $gene2tx_ref->{$gene}->{$tx} = 1;
    }

# Example of input line for mapping genes onto genome:
#
# II	gene	gene	9197458	9207148	.	+	.	Gene "WBGene00002299" ; Position "1.09195" ; Locus "let-23"
# 
# At least one case of a gene where this won't work and I need the data from a different line:
# IV	ncRNA	ncRNA_primary_transcript	8428601	8428726	.	-	.	Transcript "T26A8.6" ; Gene "WBGene00195245" ;

    # Not elsif, because that precludes using the same input line for both $tx-$gene links and $chr/$ori data!
    # Which, for at least one pesky ncRNA gene, proved necessary...
    if ( $input =~ / \A 
                        (\S+) 
                        (?: \t gene \t gene | \t ncRNA \t ncRNA_primary_transcript ) 
                        \t (\d+) 
                        \t (\d+) 
                        \t [^\t]* 
                        \t ([+]|[-]) 
                        \t [^\t]* 
                        \t [^\t]* \s* Gene \s+ \" ([^\"\s]+) \"
                     /xms ) { 
        $chr  = $1;
        $nt1  = $2;
        $nt2  = $3;
        $ori  = $4;
        $gene = $5;
                        
        # Enforce correct nt order:
        if ( $nt1 >= $nt2 ) { 
            die "Incorrect coordinate order in: $input\n";
        }

        # Work only on those chromosomes belonging to a specified OK list:
        if ( (@chromosomes) and (! exists $OK_chr{$chr} ) ) {
            next LOOP1;
        }

        # Enforce familiarity and known length of chromosome/contig:
        if (! exists $chrs{$chr} ) {
            die "Chromosome/contig $chr, cited in this input line:\n",
                "$input\n",
                "was not previously recorded in a Sequence length line in the GFF.\n",
                ;
        }

        # Record coords. and ori.
        $info_ref->{'gene'}->{$gene}->{'chr'}   = $chr;
        $info_ref->{'gene'}->{$gene}->{'start'} = $nt1;
        $info_ref->{'gene'}->{$gene}->{'end'}   = $nt2;
        $info_ref->{'gene'}->{$gene}->{'ori'}   = $ori;

        # Map gene spans onto genome:
        foreach my $nt ($nt1..$nt2) { 
            $info_ref->{'genome'}->{$chr}->{'coord'}->{$nt}->{'gene'}->{$gene} = 1;
        }
    }
}
close $GENOME or die "Can't close filehandle to genome GFF file $opts{'genome_GFF'}: $!";

# In fourth round of reading, identify starts of transcripts; link them, along with exons and introns, to genes; and map all onto the genome.
open $GENOME, '<', "$opts{'genome_GFF'}" or die "Can't open genome GFF file $opts{'genome_GFF'}: $!";

LOOP2:
while (my $input = <$GENOME>) {
    chomp $input;

    # Initialize for each text line to be evaluated:
    $chr    = q{};
    $ex_int = q{};
    $nt1    = q{};
    $nt2    = q{};
    $ori    = q{};
    $tx     = q{};
    $gene   = q{};

# Example of input line, for linking transcripts to genes, and mapping their starts to the genome:
# 
# II      Coding_transcript       protein_coding_primary_transcript       9198881 9206878 
# 	.       +       .       Transcript "ZK1067.1b"
# X	curated_miRNA	miRNA_primary_transcript	14744104	14744202	
# 	.	-	.	Transcript "C05G5.6" ; Gene "WBGene00002285" ;

    # Map starts of transcripts.
    if ( $input =~ / \A 
                     (\S+) 
                     \t [^\t]* 
                     \t \S* transcript \S* 
                     \t (\d+)
                     \t (\d+) 
                     \t [^\t]*
                     \t ([+]|[-]) 
                     \t [^\t]* 
                     \t Transcript \s+ \" ([^\"\s]+) \" 
                     [^\t]* \z /xms ) { 
        $chr = $1;
        $nt1 = $2;
        $nt2 = $3;
        $ori = $4;
        $tx  = $5;

        # If the $tx belongs to a pseudogene, ban it from all further consideration:
        if ( exists $pseudo{$tx} ) {
            next LOOP2;
        }

        # Skip unwanted chromosomes:
        if ( (@chromosomes) and (! exists $OK_chr{$chr} ) ) {
            next LOOP2;
        }

        # Again, censor the second suffix of names like "Y74C9A.2.3", and "F53G12.5a.1":
        if ( $tx =~ /\A ([^\.\s]+ \. \d+ [a-z]*) \. \d+ /xms ) {
            $tx = $1;
        }

        # Get gene while enforcing sanity:
        $gene = sanity_checked($tx, $chr, $ori, $input);

        # At last, compute start of transcript, map to genome, and link to gene.
        $tx_start_nt = q{};
        if ( $ori eq '+' ) { 
            $tx_start_nt = $nt1;
        }
        if ( $ori eq '-' ) { 
            $tx_start_nt = $nt2;
        }

        # I can't, efficiently, exclude non-adjacent genes -- I need the program to tell me what they are!
        #    but I can definitely exclude genes on a chromsome with no elements at all...
        if ( exists $info_ref->{'in_an_element'}->{'chr'}->{$chr} ) { 
            $info_ref->{'genome'}->{$chr}->{'coord'}->{$tx_start_nt}->{'tx_start_nt'}->{$gene} = 1;
        }
    }

# Example of input lines for linking introns or exons to genes and mapping them to the genome:
# 
# II	curated	intron	9199288	9199372	.	+	
# 	.	CDS "ZK1067.1b" ; Confirmed_cDNA X57767 ; Confirmed_EST RST5_373116
# II	curated	exon	9199129	9199287	.	+	.	CDS "ZK1067.1b"
# 
# I	Coding_transcript	exon	10413	10585	.	+	.	Transcript "Y74C9A.2.3"

    elsif ( $input =~ / \A
                        (\S+)
                        \t [^\t]*
                        \t \S* (exon|intron) \S*
                        \t (\d+)
                        \t (\d+)
                        \t [^\t]*
                        \t ([+]|[-])
                        \t [^\t]*
                        \t (?:CDS|Transcript) \s+ \" ([^\"\s]+) \"
                        [^\t]* \z /xms ) {
        $chr    = $1;
        $ex_int = $2;
        $nt1    = $3;
        $nt2    = $4;
        $ori    = $5;
        $tx     = $6;

        # If the $tx belongs to a pseudogene, ban it from all further consideration:
        if ( exists $pseudo{$tx} ) {
            next LOOP2;
        }

        # Skip unwanted chromosomes:
        if ( (@chromosomes) and (! exists $OK_chr{$chr} ) ) {
            next LOOP2;
        }

        # For transcripts, do any needed name truncations, and then...
        if ( $tx =~ /\A ([^\.\s]+ \. \d+ [a-z]*) \. \d+ /xms ) {
            $tx = $1;   
        }

        # Get gene while enforcing sanity:
        $gene = sanity_checked($tx, $chr, $ori, $input);

        # Map exon and intron spans onto genome, and link them to relevant gene(s) :
        foreach my $nt ($nt1..$nt2) {
            if ( exists $info_ref->{'in_an_element'}->{'chr'}->{$chr}->{'nt'}->{$nt} ) { 
                $info_ref->{'genome'}->{$chr}->{'coord'}->{$nt}->{$ex_int}->{$gene} = 1;
            }
        } 
    }
}
close $GENOME or die "Can't close filehandle to genome GFF file $opts{'genome_GFF'}: $!";

# For each nt in the genome, identify the nearest 5'-ward and 3'-ward protein-coding gene:
foreach my $chr1 (sort keys %chrs) { 
    if ( (! @chromosomes) or ( exists $OK_chr{$chr1} ) ) { 
        $chr_len = $info_ref->{'genome'}->{$chr1}->{'length'};

        # First, find nearest 5'-ward transcription start site, starting from 5'-ward telomere:
        my @chr_coords = (1..$chr_len);
        my @nearest_genes = ( '[five_primeward telomere]' );
        annotate_neighboring_genes( $chr1, \@chr_coords, \@nearest_genes, 'prot_coding', 'five_primeward');

        # Then, reverse the chrm. coords., and do the same thing for the 3'-ward TSSes:
        @chr_coords = reverse @chr_coords;
        @nearest_genes = ( '[three_primeward telomere]' );
        annotate_neighboring_genes( $chr1, \@chr_coords, \@nearest_genes, 'prot_coding', 'three_primeward');
    }
}

foreach my $element1 ( sort { namevalue($a) <=> namevalue($b) } keys %{ $info_ref->{'element'} } ) { 
    $chr = $info_ref->{'element'}->{$element1}->{'chr'};
    $nt1 = $info_ref->{'element'}->{$element1}->{'start'};
    $nt2 = $info_ref->{'element'}->{$element1}->{'end'};

    my %genes2exint  = ();
    my %ele2annots   = ();
    my %indiv_annots = ();
    my $gene_annot   = q{};

    # Define a single central $nt from which to describe the nearest 5'-ward and 3'-ward TSSs.
    my $center_nt = int(($nt1 + $nt2)/2); 

    # Orientations are 'five_primeward' or 'three_primeward'.

    # Get the gene identity (or, rarely, identities) and distance of nearest 5'-ward and 3'-ward transcription start sites:
    if (    (! exists $info_ref->{'genome'}->{$chr}->{'coord'}->{$center_nt}->{'direction'}->{'five_primeward'}->{'distance'} )
         or (! exists $info_ref->{'genome'}->{$chr}->{'coord'}->{$center_nt}->{'direction'}->{'five_primeward'}->{'gene'} ) ) {
        die "Failure to annotate nearest 5'-ward TSS at chr. $chr, nt $center_nt!\n";
    }
    if (    (! exists $info_ref->{'genome'}->{$chr}->{'coord'}->{$center_nt}->{'direction'}->{'three_primeward'}->{'distance'} )
         or (! exists $info_ref->{'genome'}->{$chr}->{'coord'}->{$center_nt}->{'direction'}->{'three_primeward'}->{'gene'} ) ) {
        die "Failure to annotate nearest 5'-ward TSS at chr. $chr, nt $center_nt!\n";
    }

    my $five_primeward_dist      = $info_ref->{'genome'}->{$chr}->{'coord'}->{$center_nt}->{'direction'}->{'five_primeward'}->{'distance'};
    my @five_primeward_genes     = sort keys %{ $info_ref->{'genome'}->{$chr}->{'coord'}->{$center_nt}->{'direction'}->{'five_primeward'}->{'gene'} };
    @five_primeward_genes        = map { $gene2name{$_} . q{ [ori} . $info_ref->{'gene'}->{$_}->{'ori'} . q{]} } @five_primeward_genes;
    my $five_primeward_genelist  = join q{/}, @five_primeward_genes;
    my $five_primeward_tss_annot = "5'-ward TSS: $five_primeward_genelist, $five_primeward_dist nt";
            
    my $three_primeward_dist      = $info_ref->{'genome'}->{$chr}->{'coord'}->{$center_nt}->{'direction'}->{'three_primeward'}->{'distance'};
    my @three_primeward_genes     = sort keys %{ $info_ref->{'genome'}->{$chr}->{'coord'}->{$center_nt}->{'direction'}->{'three_primeward'}->{'gene'} };
    @three_primeward_genes        = map { $gene2name{$_} . q{ [ori} . $info_ref->{'gene'}->{$_}->{'ori'} . q{]} } @three_primeward_genes;
    my $three_primeward_genelist  = join q{/}, @three_primeward_genes;
    my $three_primeward_tss_annot = "3'-ward TSS: $three_primeward_genelist, $three_primeward_dist nt";

    # Get a list of annotations for *each* $nt of element; condense; note "Discordant" if there's more than one.
    foreach my $nt_ele ($nt1..$nt2) { 
        my @nt_ele_g_exints = ();
        if ( exists $info_ref->{'genome'}->{$chr}->{'coord'}->{$nt_ele}->{'gene'} ) { 
            # For each gene of interest, get the exint status.
            my @wbgenes1 = sort keys %{ $info_ref->{'genome'}->{$chr}->{'coord'}->{$nt_ele}->{'gene'} };
            foreach my $nt_ele_gene (@wbgenes1) {
                @nt_ele_g_exints = ();
                if ( exists $info_ref->{'genome'}->{$chr}->{'coord'}->{$nt_ele}->{'exon'}->{$nt_ele_gene} ) {
                    push @nt_ele_g_exints, 'exon';
                }
                if ( exists $info_ref->{'genome'}->{$chr}->{'coord'}->{$nt_ele}->{'intron'}->{$nt_ele_gene} ) {
                    push @nt_ele_g_exints, 'intron';
                }
                my $nt_ele_g_exint_text = join '/', @nt_ele_g_exints;
                $genes2exint{$nt_ele_gene} = $nt_ele_g_exint_text;
            }

            # Generate a human-named list of genes with their orientations and exint status.
            my @genes1 = map { $gene2name{$_}
                               . " ($gene_status{$_}) " 
                               . '[ori'
                               . $info_ref->{'gene'}->{$_}->{'ori'} 
                               . q{, }
                               . $genes2exint{$_}
                               . q{]} 
                             }
                             @wbgenes1;
            foreach my $indiv_annot1 (@genes1) { 
                $indiv_annots{$indiv_annot1} = 1;
            }
            $gene_annot = join '; ', @genes1;
            $ele2annots{$gene_annot} = 1;
        }
    }

    my $annot_count = keys %ele2annots;
    my @final_annot_list = ();
    my $final_gene_annot = q{};
    if ( $annot_count >= 1 ) { 
        @final_annot_list = sort keys %indiv_annots;
        $final_gene_annot = join '; ', @final_annot_list;
        if ( $annot_count >= 2 ) {
            $final_gene_annot = 'Discordant: ' . $final_gene_annot;
        }
    }

    my $ele_coords =   $info_ref->{'element'}->{$element1}->{'chr'} 
                     . q{:} 
                     . $info_ref->{'element'}->{$element1}->{'start'} 
                     . q{..} 
                     . $info_ref->{'element'}->{$element1}->{'end'}
                     ;

    print $element1, 
          "\t", 
          $ele_coords,
          "\t",
          $info_ref->{'element'}->{$element1}->{'seq_length'},
          "\t",
          $info_ref->{'element'}->{$element1}->{'seq_CG'}, 
          "\t",
          $info_ref->{'element'}->{$element1}->{'best_motif_desc'},
          "\t",
          $final_gene_annot, 
          "\t", 
          $five_primeward_tss_annot, 
          "\t", 
          $three_primeward_tss_annot, 
          "\n", 
          ;
}

sub sanity_checked { 
    my ($_tx, $_chr, $_ori, $_input) = @_;

    if (! exists $chrs{$_chr} ) {
        die "Chromosome/contig $_chr, cited in this input line:\n",
            "$_input\n",
            "was not previously recorded in a Sequence length line in the GFF.\n",
            ;
    }
    if (! exists $tx2gene{$_tx}) {
        die "Failed to previously link any gene for the transcript $_tx!\n";
    }

    my $_gene = $tx2gene{$_tx};

    if (! exists $info_ref->{'gene'}->{$_gene}->{'chr'} ) {
        die "Failed to previously record any chromosome for gene $_gene, transcript $_tx on input line: $_input\n";
    }
    if (! exists $info_ref->{'gene'}->{$_gene}->{'ori'} ) {
        die "Failed to previously record any orientation for gene $_gene, transcript $_tx on input line: $_input\n";
    }
    if (    ( $_chr ne $info_ref->{'gene'}->{$_gene}->{'chr'} )
          or ( $_ori ne $info_ref->{'gene'}->{$_gene}->{'ori'} ) ) {
        die "Transcript $_tx, found in this input line:\n",
            "$_input\n",
            "has chromosome $_chr and orientation $_ori, but its parent gene ($_gene)\n",
            "instead has chromosome $info_ref->{'gene'}->{$_gene}->{'chr'} and orientation $info_ref->{'gene'}->{$_gene}->{'ori'}.\n",
            ;
    }
    return $_gene;
}

sub annotate_neighboring_genes { 
    my ( $_chr, $_chr_coords_ref, $_nearest_genes_ref, $_wanted_gene_type, $_direction) = @_;
    my $_dist_to_tx_start = 0;
    my @_nearest_genes = @{ $_nearest_genes_ref };
    foreach my $_nt ( @{ $_chr_coords_ref } ) {
        if (! exists $info_ref->{'genome'}->{$_chr}->{'coord'}->{$_nt}->{'tx_start_nt'} ) {
            $_dist_to_tx_start++;
        }
        else {
            my @_possible_nearest_genes 
                = map { accept_gene_type($_, $_wanted_gene_type) } 
                  sort keys %{ $info_ref->{'genome'}->{$_chr}->{'coord'}->{$_nt}->{'tx_start_nt'} };

            if (@_possible_nearest_genes) {
                $_dist_to_tx_start = 0;
                @_nearest_genes = @_possible_nearest_genes;
            }
        }
        # Don't actually bother to *record* values that will never actually be used!
        if ( exists $info_ref->{'in_an_element'}->{'chr'}->{$_chr}->{'nt'}->{$_nt} ) { 
            $info_ref->{'genome'}->{$_chr}->{'coord'}->{$_nt}->{'direction'}->{$_direction}->{'distance'} = $_dist_to_tx_start;
            foreach my $_neighboring_gene (@_nearest_genes) {
                $info_ref->{'genome'}->{$_chr}->{'coord'}->{$_nt}->{'direction'}->{$_direction}->{'gene'}->{$_neighboring_gene} = 1;
            }
        }
    }
    # Pro forma:
    return;
}               

sub accept_gene_type { 
    my ($_gene_id, $_gene_type) = @_;
    if (! exists $gene_status{$_gene_id} ) { 
        die "Can't decide what class of gene $_gene_id belongs to, thus can't decide whether to pass it on or not in sub accept_gene_type!\n";
    }
    elsif ( $_gene_type !~ /\S/xms ) { 
        die "Unspecified gene status for gene $_gene_id in sub accept_gene_type!\n";
    }
    elsif ( $gene_status{$_gene_id} eq $_gene_type ) { 
        return $_gene_id;
    }
    else { 
        return;
    }
}

sub namevalue { 
    # Note that this presupposes that no $_name1_2 will ever be bigger than 999_999; if *that* assumption fails [!], so does the script.
    my $_name1   = $_[0];
    my $_name1_1 = q{};
    my $_name1_2 = q{};
    my $_sortvalue = 0;
    if ( $_name1 !~ /\A (?: \d+\- | chr(?:I|II|III|IV|V|X)\. ) \d+ \z /xms ) {
        die "Can't parse name $_name1 in sorting!\n";
    }
    if ($_name1 =~ /\A ( \d+\- | chr(?:I|II|III|IV|V|X)\. )  (\d+) \z /xms ) {
        $_name1_1 = $1;
        $_name1_2 = $2;
    }
    if ( $_name1_1 =~ /\A chr(I|II|III|IV|V|X) \. \z /xms ) { 
        $_name1_1 = $1;
        $_name1_1 = arabic($_name1_1);
    }
    if ( $_name1_1 =~ /\A (\d+) \- \z /xms ) { 
        $_name1_1 = $1;
    }
    if ( $_name1_2 > 999_999 ) { 
        die "Can't compare name $_name1 in motif-sorting, because its second term ($_name1_2) is greater than 999,999!\n";
    }
    $_name1_1 = $_name1_1 * 1_000_000;
    $_sortvalue = $_name1_1 + $_name1_2;
    return $_sortvalue;
}

sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}

