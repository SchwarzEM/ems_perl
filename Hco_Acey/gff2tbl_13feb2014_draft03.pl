#!/usr/bin/env perl

# gff2tbl.pl -- Erich Schwarz <ems@emstech.org>, 2/12/2014.
# Purpose: given a GFF-like file, export a gene feature table file acceptable to NCBI Sequin; order the records in the same order as the associated FASTA.
# Note: this version has been coded to only accept AUGUSTUS outputs in GTF2 format; later versions should be expanded to handle classic GFF2 or GFF3.

use strict;
use warnings;
use Getopt::Long;

my $data_ref;

my @gff_files = ();
my $fasta     = q{};
my $locus     = q{};
my $gff_type  = q{};
my $help;

my %ok_type = ( 
    augustus => 'aug',
    aug      => 'aug',
    gtf2     => 'aug',
);

my @source_seqs = ();

GetOptions ( 'gff_files=s{,}' => \@gff_files,
             'fasta=s'        => \$fasta,
             'locus=s'        => \$locus,
             'type=s'         => \$gff_type,
             'help'           => \$help,   );

if ( $help or (! @gff_files) or (! $fasta ) or (! $locus ) or (! $gff_type ) or (! exists $ok_type{$gff_type} ) ) { 
    die "Format: gff2tbl.pl\n",
        "    --gff|-g      <input stream/files in GFF-like format>\n",
        "    --fasta|-f    [FASTA file: required to verify each sequence named in GFF(s), and to provide order for output]\n",
        "    --locus|-l    [gene-to-locus table, with variable and arbitrary mappings provided by the user]\n",
        "    --type|-t     [GFF type: currently, only 'augustus|aug' or 'gtf2'; planned additions are 'gff2' and 'gff3']\n",
        "    --help|-h     [print this message]\n",
        ;
}

open my $FASTA, '<', $fasta or die "Can't open FASTA file $fasta: $!";
while (my $input = <$FASTA>) { 
    chomp $input;
    if ( $input =~ /\A > (\S+) /xms ) { 
        my $source_seq = $1;
        if ( exists $data_ref->{'seen_in_fasta'}->{$source_seq} ) {
            die "In FASTA file $fasta, redundant sequence name $source_seq\n";
        }
        else { 
            $data_ref->{'seen_in_fasta'}->{$source_seq} = 1;
            push @{ $data_ref->{'source_seqs'} }, $source_seq;
        }
    }
    elsif ( $input =~ /\A > /xms ) { 
        die "Can't parse header in $fasta: $input\n";
    }
}
close $FASTA or die "Can't close filehandle to FASTA file $fasta: $!";

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
            # Sample input:
            # Acey_s0001_scaf AUGUSTUS        gene    248866  260706  0.88    +       .       Acey_s0001.g20
            # Use this sort of line to get *gene* nt coordinates on a source DNA sequence.
            if ( $input =~ /\A (\S+) \t [^\t]* \t gene \t (\d+) \t (\d+) \t [^\t]* \t ([+]|[-]) \t \. \t (\S+) \s* \z/xms ) { 
                my $source_seq = $1;
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
            # Sample input:
            # Acey_s0001_scaf AUGUSTUS        transcript      248866  260706  0.2     +       .       Acey_s0001.g20.t1
            # Use this sort of line to get *transcript* nt coordinates on a source DNA sequence.
            elsif ( $input =~ /\A (\S+) \t [^\t]* \t transcript \t (\d+) \t (\d+) \t [^\t]* \t ([+]|[-]) \t \. \t (\S+) \s* \z/xms ) {
                my $source_seq  = $1;
                my $start_nt    = $2;
                my $end_nt      = $3;
                my $orientation = $4;
                my $tx_id       = $5;
                $data_ref->{'source_seq'}->{$source_seq}->{'tx_id'}->{$tx_id}->{'coords'} = $start_nt . q{-} . $end_nt. q{ [} . $orientation . q{]};
            }
            # Sample input for a given transcript:
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
            # Map all of these to a transcript, and enforce a *single* gene to which each transcript maps.
            elsif ( $input =~ /\A (\S+)           # $source_seq
                                  \t [^\t]* \t 
                                  (\S+)           # $feature
                                  \t 
                                  (\d+)           # $start_nt
                                  \t 
                                  (\d+)           # $end_nt
                                  \t [^\t]* \t 
                                  ([+]|[-])       # $orientation
                                  \t 
                                  (\S+)           # $phase
                                  \t 
                                  transcript_id [ ] \" ([^\s\"]+) \" ; [ ] gene_id [ ] \" ([^\s\"]+) \" ;  # $tx_id, then $gene_id
                                  \s* 
                               \z/xms ) {
                 my $source_seq  = $1;
                 my $feature     = $2;
                 my $start_nt    = $3;
                 my $end_nt      = $4;
                 my $ori         = $5;   # orientation
                 my $phase       = $6;
                 my $tx_id       = $7;
                 my $gene_id     = $8;
                 if (     ( exists $data_ref->{'source_seq'}->{$source_seq}->{'tx_id'}->{$tx_id}->{'gene_id'}      ) 
                      and ( $gene_id ne $data_ref->{'source_seq'}->{$source_seq}->{'tx_id'}->{$tx_id}->{'gene_id'} ) ) {
                     die "In source sequence $source_seq, the transcript $tx_id is being mapped",
                         " to $gene_id despite a prior mapping to",
                         " $data_ref->{'source_seq'}->{$source_seq}->{'tx_id'}->{$tx_id}->{'gene_id'}, in text: $input\n",
                         ;
                 }
                 $data_ref->{'source_seq'}->{$source_seq}->{'tx_id'}->{$tx_id}->{'gene_id'} = $gene_id;
                 $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'tx_id'}->{$tx_id}->{'feature'}->{$feature}->{'start_nt'}->{$start_nt}->{'end_nt'} = $end_nt;
                 $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'tx_id'}->{$tx_id}->{'feature'}->{$feature}->{'start_nt'}->{$start_nt}->{'ori'} = $ori;
                 $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'tx_id'}->{$tx_id}->{'feature'}->{$feature}->{'start_nt'}->{$start_nt}->{'phase'} = $phase;
             }
             else { 
                 die "Can't parse input line from GFF-like file $gff_file (with format $gff_type): $input\n";
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
        if (! exists $data_ref->{'gene_seen_in_gff'}->{$gene_id} ) { 
            die "Gene-to-locus mapping table $locus tries to map an unattested gene: $input\n";
        }
        if ( exists $data_ref->{'gene_id'}->{$gene_id}->{'locus_id'} ) {
            die "Gene-to-locus mapping table $locus tries to map gene to locus ID twice;",
                " earlier mapping is to $data_ref->{'gene_id'}->{$gene_id}->{'locus_id'},",
                " later mapping is to $locus_id\n",
                ;
        }
        $data_ref->{'gene_id'}->{$gene_id}->{'locus_id'} = $locus_id;
    }
    else { 
        die "Can't parse text from gene-to-locus mapping table $locus: $input\n";
    }
}
close $LOCUS or die "Can't close filehandle gene-to-locus mapping table $locus: $!";

@source_seqs = grep { exists $data_ref->{'source_seq'}->{$_}->{'gene_id'} } @{ $data_ref->{'source_seqs'} };

foreach my $source_seq (@source_seqs) {
    print ">Features $source_seq Table_1\n";
    my @genes = sort keys %{ $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'} };
    foreach my $gene_id (@genes) { 
        if (! exists $data_ref->{'gene_id'}->{$gene_id}->{'locus_id'} ) {
            die "Can't find locus ID for gene $gene_id\n";
        }
        my $gene_start_nt = $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'start_nt'};
        my $gene_end_nt   = $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'end_nt'};
        my $locus_id      = $data_ref->{'gene_id'}->{$gene_id}->{'locus_id'};
        # Note: NCBI has a funny idea that *every single gene* should be listed as 'partial' if the mRNA has not been fully mapped.
        # Every single one of our gene predictions starts and stops at coding residues (start and stop codons).  So, they are absolutely, all of them, 'partial'.
        print q{<}, $gene_start_nt, "\t", q{>}, $gene_end_nt, "\tgene\n";
        print "\t\t\tgene\t$gene_id\n";
        print "\t\t\tlocus_tag\t$locus_id\n";
        my @txs = sort keys %{ $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'tx_id'} };
        my $tx_count = @txs;
        my $note_on_alt_splicing = q{};
        if ($tx_count >= 2 ) {
            $note_on_alt_splicing = 'note alternatively spliced';
        }
        foreach my $tx_id (@txs) { 
            my $product_desc   = "hypothetical protein";
            my $ncbi_id_prefix = 'gnl|ncbi|';

            my $mrna_id        = $ncbi_id_prefix . $tx_id . '.mrna';

            my $protein_id     = $tx_id;
            $protein_id        =~ tr/[a-z]/[A-Z]/;
            $protein_id        = $ncbi_id_prefix . $protein_id . '.prot';

            my @features = grep { $_ ne 'CDS' } sort keys %{ $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'tx_id'}->{$tx_id}->{'feature'} };
            my @start_cds_nts = sort { $a <=> $b }
                                keys %{ $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'tx_id'}->{$tx_id}->{'feature'}->{'CDS'}->{'start_nt'} };
            my @exon_coords = ();  
            foreach my $start_cds_nt (@start_cds_nts) { 
                my $end_cds_nt = $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'tx_id'}->{$tx_id}->{'feature'}->{'CDS'}->{'start_nt'}->{$start_cds_nt}->{'end_nt'};
                push @exon_coords, "$start_cds_nt\t$end_cds_nt";
            }
            my $exon_coord_text = join "\n", @exon_coords;
            my $mrna_coord_text = $exon_coord_text;
            $mrna_coord_text    = '>' . $mrna_coord_text;
            if ( $mrna_coord_text =~ /\A (.+\t) (\S+) \z/xms ) { 
                my $most_text = $1;
                my $last_text = $2;
                $mrna_coord_text = $most_text . q{>} . $last_text;
            }

            print "mRNA\n";
            print "$mrna_coord_text\n";
            print "\t\t\tproduct\t$product_desc\n";
            print "\t\t\tprotein_id\t$protein_id\n";
            print "\t\t\ttranscript_id\t$mrna_id\n";


            my $codon_start_text = q{};
            if ( ! exists $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'tx_id'}->{$tx_id}->{'feature'}->{'start_codon'} ) {
                $codon_start_text = 'X';
            }

            my $gene_ori = $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'ori'};
            if (    
                 ( ( ! exists $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'tx_id'}->{$tx_id}->{'feature'}->{'start_codon'} ) and ( $gene_ori eq '+' ) ) 
                 or 
                 ( ( ! exists $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'tx_id'}->{$tx_id}->{'feature'}->{'stop_codon'} ) and ( $gene_ori eq '-' ) ) 
               ) {
                $exon_coord_text = q{>} . $exon_coord_text;
            } 
            if (    
                 ( ( ! exists $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'tx_id'}->{$tx_id}->{'feature'}->{'start_codon'} ) and ( $gene_ori eq '-' ) ) 
                 or 
                 ( ( ! exists $data_ref->{'source_seq'}->{$source_seq}->{'gene_id'}->{$gene_id}->{'tx_id'}->{$tx_id}->{'feature'}->{'stop_codon'} ) and ( $gene_ori eq '+' ) ) 
               ) {
                if ( $exon_coord_text =~ /\A (.+\t) (\S+) \z/xms ) {
                    my $most_text = $1;
                    my $last_text = $2;
                    $exon_coord_text = $most_text . q{>} . $last_text;
                }
            }

            print "CDS\n";
            print "$exon_coord_text\n";
            print "\t\t\tproduct\t$product_desc\n";
            print "\t\t\tcodon_start\t$codon_start_text\n" if $codon_start_text;
            print "\t\t\tprotein_id\t$protein_id\n";
            print "\t\t\ttranscript_id\t$mrna_id\n";
        }
    }
}

