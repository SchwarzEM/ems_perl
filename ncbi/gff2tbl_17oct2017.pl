#!/usr/bin/env perl

# gff2tbl_17oct2017.pl -- Erich Schwarz <ems@emstech.org>, 10/11/2017.
# Purpose: given GFF3-like AUGUSTUS file, export gene feature table file acceptable to NCBI Sequin and curators in late 2017; 
#     order the records in the same order as the associated FASTA.
# Note: this version has been coded to only accept AUGUSTUS outputs in GTF2 or GFF3 format;
#     later versions should be expanded to handle generic GFF2 or GFF3.

use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils qw(natatime);

my $data_ref;

my @gff_files      = ();
my $fasta          = q{};
my $locus          = q{};
my $annot_gene     = q{};
my $annot_g2prot   = q{};
my $gff_type       = q{};
my $unique_lab_tag = q{};   # e.g. 'SchwarzCornell'
my $help;

my $fasta_seq = q{};

my %ok_type = ( 
    augustus => 'aug',
    aug      => 'aug',
    gtf2     => 'aug',
    gff3     => 'gff3',
);

my %ori2mult = (
    '+' =>  1,
    '-' => -1,
);

my %ori2start_ineq = (
    '+' => '<',
    '-' => '>',
);

my %ori2end_ineq = (
    '+' => '>',
    '-' => '<',
);

my @source_seqs = ();

GetOptions ( 'gff_files=s{,}' => \@gff_files,
             'fasta=s'        => \$fasta,
             'locus=s'        => \$locus,
             'a_g=s'          => \$annot_gene,
             'a_g2p=s'        => \$annot_g2prot,
             'type=s'         => \$gff_type,
             'uniq_lab=s'     => \$unique_lab_tag,
             'help'           => \$help,   );

if ( $help 
     or (! @gff_files) 
     or (! $fasta ) 
     or (! $locus ) 
     or (! $gff_type ) 
     or (! $unique_lab_tag ) 
     or (! exists $ok_type{$gff_type} ) 
   ) { 
    die "Format: gff2tbl_11oct2017.pl\n",
        "    --gff|-g       <input stream/files in GFF-like format>\n",
        "    --fasta|-f     [FASTA file: required to verify each sequence named in GFF(s), and to provide order for output]\n",
        "    --locus|-l     [gene-to-locus table, with variable and arbitrary mappings provided by the user]\n",
        "    --a_g          [optional: gene-to-gene-annotation table; useful for providing orthology-based aliases to genes, such as 'Acey-lin-3']\n",
        "    --a_g2p        [optional: gene-to-CDS/mRNA-annotation table; useful for adding information to all of a gene's proteins]\n",
        "    --type|-t      [GFF type: currently, only 'augustus|aug', (~equivalently)'gtf2', or 'gff3'; planned additions are 'gff2' and 'gff3']\n",
        "    --uniq_lab|-u  [unique laboratory tag, required by GenBank as filler for gnl|xxx|(etc.) ID texts]\n",
        "    --help|-h      [print this message]\n",
        ;
}

open my $FASTA, '<', $fasta or die "Can't open FASTA file $fasta: $!";
while (my $input = <$FASTA>) { 
    chomp $input;
    if ( $input =~ /\A > (\S+) /xms ) { 
        $fasta_seq = $1;
        if ( exists $data_ref->{'seen_in_fasta'}->{$fasta_seq} ) {
            die "In FASTA file $fasta, redundant sequence name $fasta_seq\n";
        }
        else { 
            $data_ref->{'seen_in_fasta'}->{$fasta_seq} = 1;
            push @{ $data_ref->{'source_seqs'} }, $fasta_seq;
        }
    }
    elsif ( $input =~ /\A > /xms ) { 
        die "Can't parse header in $fasta: $input\n";
    }
    else {
        $input =~ s/\s//;
        $data_ref->{'seq'}->{$fasta_seq}->{'residues'} .= $input;
    }
}
close $FASTA or die "Can't close filehandle to FASTA file $fasta: $!";

# For each sequence, get its residue count, so that later on, we can avoid the glitch of genes going right off the sequence ends.
my @seqs = sort keys %{ $data_ref->{'seq'} };
foreach my $seq (@seqs) { 
    my $residues = $data_ref->{'seq'}->{$seq}->{'residues'};
    my $length = length($residues);
    $data_ref->{'seq'}->{$seq}->{'length'} = $length;
}

foreach my $gff_file (@gff_files) { 
    my $GFF_FILE;
    if ($gff_file eq '-') {
        # Special case: get the stdin handle
        $GFF_FILE = *STDIN{IO};
    }
    else {
        # Standard case: open the file
        open $GFF_FILE, '<', $gff_file or die "Can't open GFF-like file $gff_file (with format $gff_type): $!\n";
    }
    while (my $input = <$GFF_FILE>) { 
        chomp $input;
        if ( $input !~ /\A [#]/xms ) { 
            if ( $input =~ /\A (\S+) \t /xms ) { 
                my $gff_seq = $1;
                if (! exists $data_ref->{'seen_in_fasta'}->{$gff_seq} ) { 
                    die "In GFF-like file $gff_file (with format $gff_type), there is a sequence name ($gff_seq)",
                        " that is absent from the FASTA file $fasta, in input line: $input\n",
                        ;
                }
            }

            # Use the following sorts of lines to get *gene* nt coordinates on a source DNA sequence.

            # Sample input for the overall span and orientation of a given gene, for AUGUSTUS GTF:
            # Acey_s0001_scaf AUGUSTUS        gene    248866  260706  0.88    +       .       Acey_s0001.g20

            my $source_seq = q{};

            if ( ( $ok_type{$gff_type} eq 'aug') 
                    and ( $input =~ /\A (\S+) \t [^\t]* \t gene \t (\d+) \t (\d+) \t [^\t]* \t ([+]|[-]) \t \. \t (\S+) \s* \z/xms ) ) {
                $source_seq    = $1;
                my $start_nt   = $2;
                my $end_nt     = $3;
                my $ori        = $4;
                my $gene_id    = $5;
                $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'start_nt'} = $start_nt;
                $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'end_nt'}   = $end_nt;
                $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'ori'}      = $ori;
                # Also do this, to provide a later error-check for the gene-to-locus table.
                $data_ref->{'gene_seen_in_gff'}->{$gene_id} = 1;
            }

            # Sample input for the overall span and orientation of a given gene, for AUGUSTUS GFF3:
            # nigoni.pc_2016.07.11_035	AUGUSTUS	gene	29497	34653	1.07	-	.	ID=nigoni.pc_2016.07.11_035.g28772
            # nigoni.pc_2016.07.11_035	AUGUSTUS	gene	42148	44346	0.36	-	.	ID=nigoni.pc_2016.07.11_035.g28773

            elsif ( ( $ok_type{$gff_type} eq 'gff3')
                    and ( $input =~ /\A (\S+) \t [^\t]* \t gene \t (\d+) \t (\d+) \t [^\t]* \t ([+]|[-]) \t \. \t ID= (\S+) \s* \z/xms ) ) {
                $source_seq    = $1;  
                my $start_nt   = $2;
                my $end_nt     = $3;
                my $ori        = $4; 
                my $gene_id    = $5;
                $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'start_nt'} = $start_nt;
                $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'end_nt'}   = $end_nt;
                $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'ori'}      = $ori;
                # Also do this, to provide a later error-check for the gene-to-locus table.
                $data_ref->{'gene_seen_in_gff'}->{$gene_id} = 1;
            }

            # Sample input for the overall span and orientation of a given transcript --
            #
            # In AUGUSTUS GTF2:
            # Acey_s0001_scaf AUGUSTUS        transcript      248866  260706  0.2     +       .       Acey_s0001.g20.t1
            # 
            # In AUGUSTUS GFF3:
            # nigoni.pc_2016.07.11_035	AUGUSTUS	transcript	31620	33210	0.35	-	.	ID=nigoni.pc_2016.07.11_035.g28772.t1;Parent=nigoni.pc_2016.07.11_035.g28772
            #
            # In GTF2, we don't need this (because CDS lines carry gene data.
            # But, in the more strictly formatted GFF3, we need this line so that we can map transcripts to genes!

            elsif ( ( $ok_type{$gff_type} eq 'gff3') 
                    and
                    ( $input =~ /\A (\S+)          # $source_seq
                                    \t [^\t]* \t
                                    transcript     # feature
                                    \t
                                    \d+            # start_nt
                                    \t
                                    \d+            # end_nt
                                    \t [^\t]* \t
                                    (?:[+]|[-])    # orientation
                                    \t
                                    \S+            # phase
                                    \t
                                    (.*)           # tx-to-gene annotation
                                 \z/xms )
                  ) {
                my $source_seq  = $1;
                my $tx_info = $2;

                my $tx_id       = q{};
                my $gene_id     = q{};

                # ID=nigoni.pc_2016.07.11_035.g28772.t1;Parent=nigoni.pc_2016.07.11_035.g28772
                if ( $tx_info =~ / ID= (\S+) ;  Parent= (\S+) /xms ) {
                    $tx_id   = $1;
                    $gene_id = $2;

                    # Enforce mapping of each transcript to a single, unique gene.
                    if (    ( exists $data_ref->{'source_seq'}->{$source_seq}->{'tx_id'}->{$tx_id}->{'gene_id'}      )
                        and ( $gene_id ne $data_ref->{'source_seq'}->{$source_seq}->{'tx_id'}->{$tx_id}->{'gene_id'} ) ) {
                        die "In source sequence $source_seq, the transcript $tx_id is being mapped",
                            " to $gene_id despite a prior mapping to",
                            " $data_ref->{'source_seq'}->{$source_seq}->{'tx_id'}->{$tx_id}->{'gene_id'}, in text: $input\n",
                            ;
                    }
                
                    # For each transcript, map it to a single gene.
                    $data_ref->{'source_seq'}->{$source_seq}->{'tx_id'}->{$tx_id}->{'gene_id'} = $gene_id;
                }
                else {
                    die "For AUGUSTUS/GFF3, cannot map transcript \"$tx_id\" to gene in: $input\n";
                }
            }

            # Sample input for the details of a given transcript --
            # 
            # For AUGUSTUS GTF2:
            # 
            # Acey_s0001_scaf AUGUSTUS        start_codon     248866  248868  .       +       0       transcript_id "Acey_s0001.g20.t1"; gene_id "Acey_s0001.g20";
            # Acey_s0001_scaf AUGUSTUS        intron  248953  249050  0.99    +       .       transcript_id "Acey_s0001.g20.t1"; gene_id "Acey_s0001.g20";
            # Acey_s0001_scaf AUGUSTUS        intron  249159  258137  1       +       .       transcript_id "Acey_s0001.g20.t1"; gene_id "Acey_s0001.g20";
            # Acey_s0001_scaf AUGUSTUS        intron  258237  260571  0.9     +       .       transcript_id "Acey_s0001.g20.t1"; gene_id "Acey_s0001.g20";
            # Acey_s0001_scaf AUGUSTUS        CDS     248866  248952  0.22    +       0       transcript_id "Acey_s0001.g20.t1"; gene_id "Acey_s0001.g20";
            # Acey_s0001_scaf AUGUSTUS        CDS     249051  249158  1       +       0       transcript_id "Acey_s0001.g20.t1"; gene_id "Acey_s0001.g20";
            # Acey_s0001_scaf AUGUSTUS        CDS     258138  258236  1       +       0       transcript_id "Acey_s0001.g20.t1"; gene_id "Acey_s0001.g20";
            # Acey_s0001_scaf AUGUSTUS        CDS     260572  260706  0.9     +       0       transcript_id "Acey_s0001.g20.t1"; gene_id "Acey_s0001.g20";
            # Acey_s0001_scaf AUGUSTUS        stop_codon      260704  260706  .       +       0       transcript_id "Acey_s0001.g20.t1"; gene_id "Acey_s0001.g20";
            # 
            # For AUGUSTUS GFF3:
            # 
            # nigoni.pc_2016.07.11_035	AUGUSTUS	stop_codon	31620	31622	.	-	0	Parent=nigoni.pc_2016.07.11_035.g28772.t1
            # nigoni.pc_2016.07.11_035	AUGUSTUS	intron	31850	31951	0.97	-	.	Parent=nigoni.pc_2016.07.11_035.g28772.t1
            # nigoni.pc_2016.07.11_035	AUGUSTUS	intron	32018	32671	0.97	-	.	Parent=nigoni.pc_2016.07.11_035.g28772.t1
            # nigoni.pc_2016.07.11_035	AUGUSTUS	intron	32747	33185	0.54	-	.	Parent=nigoni.pc_2016.07.11_035.g28772.t1
            # nigoni.pc_2016.07.11_035	AUGUSTUS	CDS	31620	31849	0.69	-	2	ID=nigoni.pc_2016.07.11_035.g28772.t1.cds;Parent=nigoni.pc_2016.07.11_035.g28772.t1
            # nigoni.pc_2016.07.11_035	AUGUSTUS	CDS	31952	32017	0.97	-	2	ID=nigoni.pc_2016.07.11_035.g28772.t1.cds;Parent=nigoni.pc_2016.07.11_035.g28772.t1
            # nigoni.pc_2016.07.11_035	AUGUSTUS	CDS	32672	32746	1	-	2	ID=nigoni.pc_2016.07.11_035.g28772.t1.cds;Parent=nigoni.pc_2016.07.11_035.g28772.t1
            # nigoni.pc_2016.07.11_035	AUGUSTUS	CDS	33186	33210	0.54	-	0	ID=nigoni.pc_2016.07.11_035.g28772.t1.cds;Parent=nigoni.pc_2016.07.11_035.g28772.t1
            # nigoni.pc_2016.07.11_035	AUGUSTUS	start_codon	33208	33210	.	-	0	Parent=nigoni.pc_2016.07.11_035.g28772.t1

            # Map all of these to a transcript; 
            #     enforce a *single* gene to which each transcript maps;
            #     identify individual CDSes by their start nts;
            #     keep track, for each CDS, its end nt, orientation, and phase.

            # Parse AUGUSTUS GTF2 or AUGUSTUS GFF3 gene body data:

            elsif ( ( ( $ok_type{$gff_type} eq 'aug') or ( $ok_type{$gff_type} eq 'gff3') ) 
                    and 
                    ( $input =~ /\A (\S+)           # $source_seq
                                    \t [^\t]* \t 
                                    (start_codon|CDS|intron|stop_codon)   # $feature
                                    \t 
                                    (\d+)           # $start_nt
                                    \t 
                                    (\d+)           # $end_nt
                                    \t [^\t]* \t 
                                    ([+]|[-])       # $orientation
                                    \t 
                                    (\S+)           # $phase
                                    \t
                                    (.*)
                                 \z/xms ) 
                  ) {
                my $source_seq  = $1;
                my $feature     = $2;
                my $start_nt    = $3;
                my $end_nt      = $4;
                my $ori         = $5;   # orientation
                my $phase       = $6;
                my $tx_info     = $7;

                my $tx_id       = q{};
                my $gene_id     = q{};

                # AUGUSTUS GTF2-specific parsing of transcript and gene ID text:
                if ( $ok_type{$gff_type} eq 'aug') {
                    if ( $tx_info =~ / transcript_id [ ] \" ([^\s\"]+) \" ; [ ] gene_id [ ] \" ([^\s\"]+) \" \s* /xms ) { 
                        $tx_id   = $1;
                        $gene_id = $2;
                    }
                    else {
                        die "For AUGUSTUS/GTF, cannot parse tx/gene info in: $input\n";
                    }

                    # Enforce mapping of each transcript to a single, unique gene.
                    if (    ( exists $data_ref->{'source_seq'}->{$source_seq}->{'tx_id'}->{$tx_id}->{'gene_id'}      )
                        and ( $gene_id ne $data_ref->{'source_seq'}->{$source_seq}->{'tx_id'}->{$tx_id}->{'gene_id'} ) ) {
                        die "In source sequence $source_seq, the transcript $tx_id is being mapped",
                            " to $gene_id despite a prior mapping to",
                            " $data_ref->{'source_seq'}->{$source_seq}->{'tx_id'}->{$tx_id}->{'gene_id'}, in text: $input\n",
                            ;
                    }

                    # For each transcript, map it to a single gene.
                    $data_ref->{'source_seq'}->{$source_seq}->{'tx_id'}->{$tx_id}->{'gene_id'} = $gene_id;
                } 

                # AUGUSTUS GFF3-specific parsing of transcript (and, implicitly, gene) ID text:
                if ( $ok_type{$gff_type} eq 'gff3') {
                    if ( $tx_info =~ / Parent= (\S+) /xms ) {
                        $tx_id = $1;
                        if ( exists $data_ref->{'source_seq'}->{$source_seq}->{'tx_id'}->{$tx_id}->{'gene_id'} ) {
                            $gene_id = $data_ref->{'source_seq'}->{$source_seq}->{'tx_id'}->{$tx_id}->{'gene_id'};
                        }
                        else { 
                            die "For AUGUSTUS/GFF3, cannot map transcript \"$tx_id\" to gene in: $input\n";
                        }
                    }
                    else {
                        die "For AUGUSTUS/GFF3, cannot parse tx/gene info in: $input\n";
                    }
                }

                # Identity each feature (e.g., CDS, exon, intron, etc.) of the transcript by its start_nt, 
                #     and track the linked end_nt, ori, and phase:
                $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'tx_id'}->{$tx_id}->{'feature'}->{$feature}->{'start_nt'}->{$start_nt}->{'end_nt'} = $end_nt;
                $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'tx_id'}->{$tx_id}->{'feature'}->{$feature}->{'start_nt'}->{$start_nt}->{'ori'}    = $ori;
                $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'tx_id'}->{$tx_id}->{'feature'}->{$feature}->{'start_nt'}->{$start_nt}->{'phase'}  = $phase;

                # For start and stop codons, we have a particular need to record the exact first nt of the start, and the exact last nt of the stop.
                if ( $feature eq 'start_codon' ) { 
                    if ( $ori eq '+' ) { 
                        $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'start_codon_nt'}->{$start_nt} = 1;
                    }
                    if ( $ori eq '-' ) {
                        $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'start_codon_nt'}->{$end_nt} = 1;
                    }
                }
                if ( $feature eq 'stop_codon' ) {
                    if ( $ori eq '+' ) {
                        $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'stop_codon_nt'}->{$end_nt} = 1;
                    }
                    if ( $ori eq '-' ) {
                        $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'stop_codon_nt'}->{$start_nt} = 1;
                    }
                }       
            }
        }
    }
    close $GFF_FILE or die "Can't close filehandle to GFF-like file $gff_file (with format $gff_type): $!";
}

open my $LOCUS, '<', $locus or die "Can't open gene-to-locus mapping table $locus: $!";
while (my $input = <$LOCUS>) {
    chomp $input;
    if ( $input =~ /\A (\S+) \t (\S+) /xms ) {
        my $gene_id  = $1;
        my $locus_id = $2;

        # Enforce existence of the gene in the GFF (and, by earlier requirement, in the FASTA as well).
        if (! exists $data_ref->{'gene_seen_in_gff'}->{$gene_id} ) { 
            die "Gene-to-locus mapping table $locus tries to map an unattested gene: $input\n";
        }

        # Enforce unique gene-to-locus mapping.
        if ( exists $data_ref->{'gene_id'}->{$gene_id}->{'locus_id'} ) {
            die "Gene-to-locus mapping table $locus tries to map gene to locus ID twice;",
                " earlier mapping is to $data_ref->{'gene_id'}->{$gene_id}->{'locus_id'},",
                " later mapping is to $locus_id\n",
                ;
        }

        # Map each gene to one locus.  Since the file is user-specified, many different mappings are possible.
        $data_ref->{'gene_id'}->{$gene_id}->{'locus_id'} = $locus_id;
    }
    # Enforce successful parsing.
    else { 
        die "Can't parse text from gene-to-locus mapping table $locus: $input\n";
    }
}
close $LOCUS or die "Can't close filehandle to gene-to-locus mapping table $locus: $!";

# Providing a table mapping genes to *gene* annotations is optional, not mandatory.
if ($annot_gene) {
    open my $ANNOT_GENE, '<', $annot_gene or die "Can't open gene-to-gene-annotation mapping table $annot_gene: $!";
    while (my $input = <$ANNOT_GENE>) {
        chomp $input;
        if ( $input =~ /\A (\S+) \t (\S+) \z/xms ) {
            my $gene_id         = $1;
            my $annotation_text = $2;

            # Note that we are no longer using the text contents of $annotation_text not merely to provide additional information, 
            #     but the *tag* for that information.
            # NCBI has started banning such things as:
            #
            # "gene_synonym\tAcey-lin-3"
            # 
            # For a list of possible gene tags, see:
            #     http://www.insdc.org/documents/feature_table.html
            #     [Scroll down to "Feature Key           gene".]

            # Enforce existence of the gene in the GFF (and, by earlier requirement, in the FASTA as well).

            if (! exists $data_ref->{'gene_seen_in_gff'}->{$gene_id} ) {
                die "Gene-to-gene-annotation mapping table $annot_gene tries to map an unattested gene: $input\n";
            }
                     
            # However, do *not* enforce unique gene-to-gene-annotation mapping; instead, we will collect all annotations and provide them as 'note X' lines later.
            # Map each gene to one or more non-redundant gene annotations.  Since the file is user-specified, many different mappings are possible.

            $data_ref->{'gene_id'}->{$gene_id}->{'gene_annotation'}->{$annotation_text} = 1;
        }
        # Enforce successful parsing.
        else {
            die "Can't parse text from gene-to-annotation mapping table $annot_gene: $input\n";
        }
    }
    close $ANNOT_GENE or die "Can't close filehandle to gene-to-annotation mapping table $annot_gene: $!";
}

# Providing a table mapping genes to their *protein* (mRNA and CDS) annotations is also optional, not mandatory.
if ($annot_g2prot) {
    open my $ANNOT_G2PROT, '<', $annot_g2prot or die "Can't open gene-to-protein-annotation mapping table $annot_g2prot: $!";
    while (my $input = <$ANNOT_G2PROT>) {
        chomp $input;
        if ( $input =~ /\A (\S+) \t (\S+\t.*) \z/xms ) {
            my $gene_id         = $1;
            my $annotation_text = $2;

            # Again, we are using the text contents of $annotation_text not merely to provide additional information, but the *tag* for that information.
            # For a list of possible gene tags, see:
            #     http://www.insdc.org/documents/feature_table.html
            #     [Scroll down to "Feature Key           CDS [or] mRNA".]

            # Enforce existence of the gene in the GFF (and, by earlier requirement, in the FASTA as well).
            if (! exists $data_ref->{'gene_seen_in_gff'}->{$gene_id} ) {
                die "Gene-to-protein-annotation mapping table $annot_g2prot tries to map an unattested gene: $input\n";
            }
                     
            # However, do *not* enforce unique gene-to-annotation mapping; instead, we will collect all annotations and provide them as 'note X' lines later.
            # Map each gene to one or more non-redundant annotations.  Since the file is user-specified, many different mappings are possible.
            $data_ref->{'gene_id'}->{$gene_id}->{'prot_annotation'}->{$annotation_text} = 1;
        }
        # Enforce successful parsing.
        else {
            die "Can't parse text from gene-to-protein-annotation mapping table $annot_g2prot: $input\n";
        }
    }
    close $ANNOT_G2PROT or die "Can't close filehandle to gene-to-protein-annotation mapping table $annot_g2prot: $!";
}

# List those sequences for which we had *some* annotation.
@source_seqs = grep { exists $data_ref->{'source_seq'}->{$_}->{'gene_id'} } @{ $data_ref->{'source_seqs'} };

# The following code generates a table, first by sequences, then by genes within it.  The order of genes is not necessarily nt-based, and could easily not be.
# If this is a concern (say, if tbl2asn stupidly chokes), it may be necessary to allow user-specific gene orders from yet another input file.

foreach my $source_seq (@source_seqs) {
    # "Feature_Table" could probably be skipped.
    print ">Features $source_seq Feature_Table\n";

    my @genes = sort keys %{ $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'} };
    foreach my $gene_id (@genes) { 
        # Require mapping to locus.
        if (! exists $data_ref->{'gene_id'}->{$gene_id}->{'locus_id'} ) {
            die "Can't find locus ID for gene $gene_id\n";
        }

        my $gene_start_nt = $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'start_nt'};
        my $gene_end_nt   = $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'end_nt'};
        my $locus_id      = $data_ref->{'gene_id'}->{$gene_id}->{'locus_id'};
        my $gene_ori      = $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'ori'};

        # If on the opposite strand, reverse the coordinates, because GenBank uses pre-modern text-based ordering of coordinates rather than '+' or '-'.
        if ( $gene_ori eq '-' ) { 
            ($gene_start_nt, $gene_end_nt) = ($gene_end_nt, $gene_start_nt);
        }

        # Note: NCBI has the idea that *every single gene* should be listed as 'partial' if the mRNA has not been fully mapped.
        # Unless otherwise noted, every single one of our gene predictions starts and stops at coding residues (start and stop codons).
        # So, they are absolutely, all of them, 'partial', according to NCBI.
        # However, in real life, this seems to produce massive failure in tbl2asn outputs; so, try giving NCBI only a little of what it wants, 
        #     when it might actually make sense.

        # This part is nice and messy; but we can reasonably assume that we have complete boundaries for a gene if there are start or 
        #     end codon nt at its furthest known coordinates.
        # However, there is one caveat: tbl2asn seems not to be able to handle a gene with '<1' or '>max_nt_of_seq' ends.
        # So, in cases where we would otherwise do that, block the application of the rule.
        my $end_seq_nt = $data_ref->{'seq'}->{$source_seq}->{'length'};

        if ( ( $gene_start_nt != 1 ) and ( $gene_start_nt != $end_seq_nt ) 
             and ( ! exists $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'start_codon_nt'}->{$gene_start_nt} ) ) { 
            $gene_start_nt = $ori2start_ineq{$gene_ori} . $gene_start_nt;
        }
        if ( ( $gene_end_nt != 1 ) and ( $gene_end_nt != $end_seq_nt )    
             and ( ! exists $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'stop_codon_nt'}->{$gene_end_nt} ) ) {
            $gene_end_nt = $ori2end_ineq{$gene_ori} . $gene_end_nt;
        }

        # NCBI has seriously changed its rules; all the gene names must now go in a single 'note' right under 'gene'.
        # Before I put text in that 'note', I need to be sure I know what the gene names I want are.

        my @ncbi_gene_names = ();
        push @ncbi_gene_names, $gene_id;
        if ( exists $data_ref->{'gene_id'}->{$gene_id}->{'gene_annotation'} ) {
            my @gene_annots = sort keys %{ $data_ref->{'gene_id'}->{$gene_id}->{'gene_annotation'} };
            foreach my $gene_annot (@gene_annots) {
                push @ncbi_gene_names, $gene_annot;
            }
        }
        my $ncbi_gene_name_note = join '; ', @ncbi_gene_names;

        # With or without added-on '<', '>', etc:
        print "$gene_start_nt\t$gene_end_nt\tgene\n";

        print "\t\t\tnote\t$ncbi_gene_name_note\n";
        print "\t\t\tlocus_tag\t$locus_id\n";

        # For the gene being considered, check to see if we will need to note alternative splicing, later on.
        my @txs = sort keys %{ $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'tx_id'} };
        my $tx_count = @txs;
        my $note_on_alt_splicing = q{};
        if ($tx_count >= 2 ) {
            $note_on_alt_splicing = 'alternatively spliced';
        }

        # Then, for each gene, work through its transcripts.
        foreach my $tx_id (@txs) { 
            my $ncbi_id_prefix = 'gnl|' . $unique_lab_tag . q{|};

            my $protein_name   = $tx_id;

            my $isoform_suffix = q{};
            if ( $protein_name =~ /\A \S+ (\.t\d+) \z/xms ) {
                $isoform_suffix = $1;
            }
            else {
                die "Can't parse protein name $protein_name\n";
            }

            my $mrna_id    = $ncbi_id_prefix . $locus_id . $isoform_suffix . '.mrna' ;

            $protein_name      =~ tr/[a-z]/[A-Z]/;

            my $protein_id = $locus_id . $isoform_suffix;
            $protein_id    =~ tr/[a-z]/[A-Z]/;
            $protein_id    = $ncbi_id_prefix . $protein_id;

            my $product_desc   = 'hypothetical protein';   # was formerly "$protein_name";

            my @start_cds_nts = sort { $a <=> $b }
                                keys %{ $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'tx_id'}->{$tx_id}->{'feature'}->{'CDS'}->{'start_nt'} };

            my @exon_coords     = ();
            my @exon_phases     = ();
            my $exon_coord_text = q{};

            foreach my $start_cds_nt (@start_cds_nts) { 
                my $end_cds_nt = $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'tx_id'}->{$tx_id}->{'feature'}->{'CDS'}->{'start_nt'}->{$start_cds_nt}->{'end_nt'};
                my $cds_phase  = $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'tx_id'}->{$tx_id}->{'feature'}->{'CDS'}->{'start_nt'}->{$start_cds_nt}->{'phase'};

                # Start by making these two lists as if the orientation was '+'; later on, we may end up reversing it.
                push @exon_coords, $start_cds_nt, $end_cds_nt;
                push @exon_phases, $cds_phase;
            }

            # Now, if the gene is antisense, the arrays get reversed:
            # Now, if the gene is antisense, the array itself (and not just its start and stop nt positions) gets reversed:
            if ( $gene_ori eq '-' ) {
                @exon_coords = reverse @exon_coords;
                @exon_phases = reverse @exon_phases;
            }

            # It makes life easier if we turn this into a line of text, with pairs of coordinates joined by '\t', and split by internal (self-splitting) '\n' characters.
            # Learning Perl, 6th. ed., p. 293:
            my $iterator = natatime 2, @exon_coords;
            while( my @pair = $iterator->() ) {
                my $exon_coord_pair = join "\t", @pair;
                $exon_coord_text .= "$exon_coord_pair\n";
            }
            chomp $exon_coord_text;

            # NCBI wants us to consider all mRNAs 'partial' here, because we are making mRNA predictions that start and stop with CDS.
            # However, I am mostly squelching this, because it appears to fail miserably in practice when fed to tbl2asn.
            #    I will only allow it for mRNAs that are truly known to be partial (because they are from genes with no observed/predicted start or stop codon).

            my $mrna_coord_text = $exon_coord_text;

            if ( ! exists $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'tx_id'}->{$tx_id}->{'feature'}->{'start_codon'} ) {
                $mrna_coord_text    = $ori2start_ineq{$gene_ori} . $mrna_coord_text;
            }

            if ( $mrna_coord_text =~ /\A (.+\t) (\S+) \z/xms ) { 
                my $most_text = $1;
                my $last_text = $2;
                if ( ! exists $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'tx_id'}->{$tx_id}->{'feature'}->{'stop_codon'} ) {
                    $mrna_coord_text = $most_text . $ori2end_ineq{$gene_ori} . $last_text;
                }
                else { 
                    $mrna_coord_text = $most_text . $last_text;
                }
            }

            # Add in 'mRNA' onto what will become the first line of $mrna_coord_text.
            if ( $mrna_coord_text =~ /\A (\S+ \t \S+) (.*) \z/xms ) {
                my $front_mrna_text = $1;
                my $back_mrna_text  = $2;
                $mrna_coord_text    = $front_mrna_text . "\tmRNA" . $back_mrna_text;
            }

            # Print mRNA text:
            print "$mrna_coord_text\n";  # The first line of this should have "\d+\t\d+\tmRNA\n".
            print "\t\t\tproduct\t$product_desc\n";
            print "\t\t\tnote\t$protein_name\n";
            print "\t\t\tnote\t$note_on_alt_splicing\n" if $note_on_alt_splicing;
            if ( exists $data_ref->{'gene_id'}->{$gene_id}->{'prot_annotation'} ) {
                my @prot_annots = sort keys %{ $data_ref->{'gene_id'}->{$gene_id}->{'prot_annotation'} };
                foreach my $prot_annot (@prot_annots) {
                    print "\t\t\t$prot_annot\n";
                }
            }
            print "\t\t\tprotein_id\t$protein_id\n";
            print "\t\t\ttranscript_id\t$mrna_id\n";

            my $codon_start_text = q{};   # Default is that this is unnecessary, ergo empty text.
            if ( ! exists $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'tx_id'}->{$tx_id}->{'feature'}->{'start_codon'} ) {
                $codon_start_text = ($exon_phases[0] + 1);
            }

            # It may actualy be a good idea to follow NCBI's silly '<' and '>' conventions: namely, when we have partial coding sequences.
            # 
            # Here, the pre-modern formatting of NCBI's tables actually helps us a bit,
            #     because start codon and stop codons are *always* at the start and stop of $exon_coord_text.

            if ( ! exists $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'tx_id'}->{$tx_id}->{'feature'}->{'start_codon'} ) { 
                $exon_coord_text = $ori2start_ineq{$gene_ori} . $exon_coord_text;
            } 
            if ( ! exists $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'tx_id'}->{$tx_id}->{'feature'}->{'stop_codon'} ) {
                if ( $exon_coord_text =~ /\A (.+\t) (\S+) \z/xms ) {
                    my $most_text = $1;
                    my $last_text = $2;
                    $exon_coord_text = $most_text . $ori2end_ineq{$gene_ori} . $last_text;
                }
            }

            # Add in 'CDS' onto what will become the first line of $exon_coord_text.   
            if ( $exon_coord_text =~ /\A (\S+ \t \S+) (.*) \z/xms ) {
                my $front_exon_text = $1;
                my $back_exon_text  = $2;
                $exon_coord_text    = $front_exon_text . "\tCDS" . $back_exon_text;
            }

            # Print CDS text:
            print "$exon_coord_text\n";  # The first line of this should have "\d+\t\d+\tCDS\n".
            print "\t\t\tproduct\t$product_desc\n";
            print "\t\t\tnote\t$protein_name\n";
            print "\t\t\tnote\t$note_on_alt_splicing\n" if $note_on_alt_splicing;
            if ( exists $data_ref->{'gene_id'}->{$gene_id}->{'prot_annotation'} ) {
                my @prot_annots = sort keys %{ $data_ref->{'gene_id'}->{$gene_id}->{'prot_annotation'} };
                foreach my $prot_annot (@prot_annots) {
                    print "\t\t\t$prot_annot\n";
                }
            }

            print "\t\t\tcodon_start\t$codon_start_text\n" if $codon_start_text;
            print "\t\t\tprotein_id\t$protein_id\n";
            print "\t\t\ttranscript_id\t$mrna_id\n";
        }
    }
}

