#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use Getopt::Long;
use File::Spec::Functions;
use Scalar::Util qw(looks_like_number);

my $genome_seq = q{};
my $sepal_set  = q{};
my @fimo_files = ();
my $print_seqs;

my $data_ref;

my $motif = q{};

my $header = "Data\t"
             . "Motif\t"
             . "All_genes\t"
             . "All_sepal_genes\t"
             . "Genomewide_hits\t"
             . "Sepalwide_hits\t"
             . "Gwide_hit_ratio\t"
             . "Sepalwide_hit_ratio\t"
             . "Sepal_enrichment"
             ;

my $help;

GetOptions ( 'genome_seq=s'    => \$genome_seq,
             'sepal_set=s'     => \$sepal_set,
             'fimo_files=s{,}' => \@fimo_files,
             'print_seqs'      => \$print_seqs,
             'help'            => \$help,        );

if ($help or (! $genome_seq ) or (! $sepal_set ) or (! @fimo_files ) ) { 
    die "Format: get_fimo_hit_rates_12mar2016.pl\n",
        "            -g|--genome_seq  [genome seq. file, from which genome-wide gene names can be extracted]\n",
        "            -s|--sepal_set   [list (not FASTA file) of sepal-wide gene names]\n",
        "            -f|--fimo_files  [1+ FIMO files, with full enough file names to identify the data type]\n",
        "            -p|--print_seqs  [optionally, print individual list, for each motif, of all the genes it detected]\n",
        "            -h|--help       [print this message]\n",
        "            [print summary table to <STDOUT>]\n",
        ;
}


open my $GEN, '<', $genome_seq;
while (my $input = <$GEN>) {
    # Note that this maps contig names like 'AT1G38630_2' to gene names like 'AT1G38630'.
    # For our purposes that is appropriate; we don't care if a gene got hit 2+ times in 2+ contigs; only if it was hit at all.
    if ( $input =~ /\A [>] (AT (?: 1|2|3|4|5|C|M) G\d+) /xms ) {
        my $seq = $1;
        $data_ref->{'genseq'}->{$seq} = 1;
    }
    elsif ( $input =~ /\A [>] /xms ) { 
        die "Cannot parse header in genome seq. $genome_seq: $input\n";
    }
}
close $GEN;

open my $SEPAL, '<', $sepal_set;
while (my $input = <$SEPAL>) {
    # Record genes only, not contig suffixes.
    if ( $input =~ /\A (AT (?: 1|2|3|4|5|C|M) G\d+) /xms ) {
        my $seq = $1;
        $data_ref->{'sepal_gene'}->{$seq} = 1;
    }
    else {
        die "Cannot parse text in sepal gene list $sepal_set: $input\n";
    }
}
close $SEPAL;

foreach my $fimo_file (@fimo_files) {
    if (! -r $fimo_file) {
        die "Cannot read putative FIMO input file: $fimo_file\n";
    }
    if ( $fimo_file =~ /\A .* \/ [^\/\s]+ \/ ([^\/\s]+) \/ fimo\.txt \z /xms ) {
        my $data_type = $1;
        open my $FIMO, '<', $fimo_file;
        $motif = q{};
        while (my $input = <$FIMO>) {
            chomp $input;
            # Once again, we take up gene names but do not distinguish contig names (if they happen to be there).
            # In this particular set of inputs, we have no p-values > 1e-05, so we do not need to impose a filter; otherwise, we would.
            # Sample input:
            # 3	AT1G38630_2	2	13	-	15.6279	4.03e-06	0.0153	AAGAAAAAGGAG
            if ( $input =~ /\A (\S+) \t (AT (?: 1|2|3|4|5|C|M) G\d+) [^\t]* \t (?: [^\t]* \t){4} \S+ \t /xms ) {
                $motif      = $1;
                my $gene    = $2;
                if (! exists $data_ref->{'genseq'}->{$gene} ) {
                    die "From FIMO file $fimo_file, cannot recognize gene \"$gene\" in: $input\n";
                }
                # We would impose a conditional based on p-values here:
                $data_ref->{'data_type'}->{$data_type}->{'motif'}->{$motif}->{'genomewide'}->{$gene} = 1;
            }
            elsif ( $input !~ /\A [#] pattern [ ] name /xms ) { 
                die "From FIMO file $fimo_file, cannot parse: $input\n";
            }
        }
        close $FIMO;

        my @motifs = sort keys %{ $data_ref->{'data_type'}->{$data_type}->{'motif'} };
        foreach my $motif_1 (@motifs) {
            my @genome_wide_hits = ();
            # Defend ourselves against cases where *zero* motif sites are found genomewide.
            if ( exists $data_ref->{'data_type'}->{$data_type}->{'motif'}->{$motif_1}->{'genomewide'} ) {
                @genome_wide_hits = sort keys %{ $data_ref->{'data_type'}->{$data_type}->{'motif'}->{$motif_1}->{'genomewide'} };
            }

            # Of the genome-wide hits, how many overlapped in any way with our sepal genes that changed expression genotypically?
            # This is a crude measurement and can miss subtle subsets, but it's also not particularly biased by any specific ideas.
            my @sepal_spec_gwide_hits = grep { ( exists $data_ref->{'sepal_gene'}->{$_} ) } @genome_wide_hits;

            my @total_genes      = sort keys %{ $data_ref->{'genseq'} };
            my @sepal_genes      = sort keys %{ $data_ref->{'sepal_gene'} };
            my @orig_genes       = sort keys %{ $data_ref->{'data_type'}->{$data_type}->{'orig_seq_set'} };

            my $genomewide_count  = @genome_wide_hits;
            my $sepalwide_count   = @sepal_spec_gwide_hits;
            my $total_gene_count  = @total_genes;
            my $total_sepal_count = @sepal_genes;

            if ( ( $total_gene_count == 0 ) or ( $total_sepal_count == 0 ) ) {
                die "Pathological situation: total gene count is $total_gene_count, total sepal gene count is $total_sepal_count\n";
            }
            my $genomewide_ratio = ($genomewide_count / $total_gene_count);
            my $sepalwide_ratio  = ($sepalwide_count  / $total_sepal_count);

            my $sepal_enrichment  = 0;
            if ( $genomewide_ratio > 0 ) {
                $sepal_enrichment = ($sepalwide_ratio / $genomewide_ratio);
            }

            # Optionally, as a possibly useful side product, print list of genomewide hits for each motif.
            # This can be given a fast check with AgriGO, and followed up if it looks promising.
            if ($print_seqs) {
                my $gwide_motif_hits = "$data_type.motif.$motif_1.genomewide_hits.txt";
                $gwide_motif_hits = safename($gwide_motif_hits);
                open my $MOT_HITS, '>', $gwide_motif_hits;
                foreach my $gwide_hit (@genome_wide_hits) {
                    print $MOT_HITS "$gwide_hit\n";
                }
                close $MOT_HITS;
            }

            # format at last possible step:
            my $print_genomewide_count = commify($genomewide_count);
            my $print_sepalwide_count  = commify($sepalwide_count);
            my $print_total_gene_count = commify($total_gene_count);
            my $print_tot_sep_g_count  = commify($total_sepal_count);

            my $print_gwide_ratio    = sprintf "%.4f", $genomewide_ratio;
            my $print_swide_ratio    = sprintf "%.4f", $sepalwide_ratio;
            my $print_s_enrichment   = sprintf "%.2f", $sepal_enrichment;

            print "$header\n" if $header;
            $header = q{};

            print "$data_type\t";
            print "$motif_1\t";
            print "$print_total_gene_count\t";
            print "$print_tot_sep_g_count\t";
            print "$print_genomewide_count\t";
            print "$print_sepalwide_count\t";
            print "$print_gwide_ratio\t";
            print "$print_swide_ratio\t";
            print "$print_s_enrichment";
            print "\n";
        }
    }
    else { 
        die "Cannot parse input: $fimo_file\n";
    }    
}

# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}

sub safename {
    my $filename = $_[0];
    my $orig_filename = $filename;
    if (-e $orig_filename) {
        my $suffix1 = 1;
        $filename = $filename . ".$suffix1";
        while (-e $filename) {
            $suffix1++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix1";
        }
    }
    return $filename;
}


