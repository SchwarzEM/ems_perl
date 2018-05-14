#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

use File::Spec::Functions;
use Scalar::Util qw(looks_like_number);

my $genome_seq = q{};
my $fimo_list  = q{};

$genome_seq = $ARGV[0] or die "Format: get_fimo_hit_rates_08mar2016.pl [genome seq. file] [fimo list file]\n";
$fimo_list  = $ARGV[1] or die "No, really -- Format: get_fimo_hit_rates_08mar2016.pl [genome seq. file] [fimo list file]\n";

my $data_ref;

my $motif = q{};

my $header = "Data\tMotif\tRedisc_founder_genes\tOrig_founder_genes\tGenomewide_hits\tAll_genes\tRedisc_ratio\tGwide_ratio\tEnrichment";

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

open my $LIST, '<', $fimo_list;
while (my $fimo_file = <$LIST> ) {
    chomp $fimo_file;
    if (! -r $fimo_file) {
        die "Cannot read putative FIMO input file: $fimo_file\n";
    }

    if ( $fimo_file =~ /\A \/mnt\/home\/emsch\/work\/2015\/adrienne\/ncDNA_motifs\/ ([^\/\s]+) \/ ( [^\/\s]+ \/ fimo\.txt ) \z /xms ) {
        my $data_type = $1;
        my $fimo_stem = $2;

        my $meme_stem = $fimo_stem;
        $meme_stem =~ s/fimo/meme/g;
        my $meme_file = catfile('/mnt/home/emsch/work/2015/adrienne/ncDNA_motifs', $data_type, $meme_stem);

        if (! -r $meme_file ) {
            die "Cannot read putative MEME input file: $meme_file\n";
        }

        open my $FIMO, '<', $fimo_file;
        $motif = q{};
        while (my $input = <$FIMO>) {
            chomp $input;
            # Once again, we take up gene names but do not distinguish contig names (if they happen to be there).
            if ( $input =~ /\A (\S+) \t (AT (?: 1|2|3|4|5|C|M) G\d+) [\t]* \t (?: [^\t]* \t){4} (\S+) \t /xms ) {
                $motif      = $1;
                my $gene    = $2;
                my $p_value = $3;

                if (! exists $data_ref->{'genseq'}->{$gene} ) {
                    die "From FIMO file $fimo_file, cannot recognize gene \"$gene\" in: $input\n";
                }
                if (! looks_like_number($p_value) ) {
                    die "From FIMO file $fimo_file, non-numerical p-value \"$p_value\" in: $input\n";
                }

                if ( $p_value <= 1e-05 ) {
                    $data_ref->{'data_type'}->{$data_type}->{'motif'}->{$motif}->{'genomewide'}->{$gene} = 1;
                }
            }
            elsif ( $input !~ /\A [#] pattern [ ] name /xms ) { 
                die "From FIMO file $fimo_file, cannot parse: $input\n";
            }
        }
        close $FIMO;

        open my $MEME, '<', $meme_file;
        # We use $motif in two ways; first, obviously, for its identity; but also, as a flag to tell us when to count gene hits!
        $motif = q{};
        while (my $input = <$MEME>) {
            chomp $input;
            if ( $input =~ / Motif [ ] (\S+) [ ] sites [ ] sorted [ ] by [ ] position [ ] p-value /xms ) {
                $motif = $1;
            }
            # Once again, we record gene names but ignore contig suffixes, *if* we have a non-zero $motif.
            elsif ( $motif and ( $input =~ /\A (AT (?: 1|2|3|4|5|C|M) G\d+) [\t]* /xms ) ) {
                my $gene = $1;
                if (! exists $data_ref->{'genseq'}->{$gene} ) {
                    die "From MEME file $meme_file, cannot record a hit for motif $motif in gene $gene\n";
                }
                elsif (! $data_ref->{'data_type'}->{$data_type}->{'motif'}->{$motif}->{'genomewide'} ) {
                    warn "Unusual and perhaps pathological situation: for data type \"$data_type\", ",
                        "motif $motif was found in original MEME search but not in later FIMO search!\n",
                        ;
                }
                else { 
                    $data_ref->{'data_type'}->{$data_type}->{'motif'}->{$motif}->{'orig_positive'}->{$gene} = 1;
                }
            }
            elsif ( $input =~ / Motif \S+ block diagrams /xms ) { 
                # Zero out the identity and also negate the "read gene hits" flag.
                $motif = q{};
            }
        }
        close $MEME;

        my @motifs = sort keys %{ $data_ref->{'data_type'}->{$data_type}->{'motif'} };
        foreach my $motif_1 (@motifs) {
            my @redisc_origs = ();
            my @orig_genes   = sort keys %{ $data_ref->{'data_type'}->{$data_type}->{'motif'}->{$motif_1}->{'orig_positive'} };

            foreach my $orig_gene (@orig_genes) {
                if ( exists $data_ref->{'data_type'}->{$data_type}->{'motif'}->{$motif_1}->{'genomewide'}->{$orig_gene} ) {
                    push @redisc_origs, $orig_gene;
                }
            }

            my @genome_wide_hits = sort keys %{ $data_ref->{'data_type'}->{$data_type}->{'motif'}->{$motif_1}->{'genomewide'} };
            my @total_genes      = sort keys %{ $data_ref->{'genseq'} };

            my $redisc_gene_count = @redisc_origs;
            my $orig_gene_count   = @orig_genes;

            my $genomewide_count  = @genome_wide_hits;
            my $total_gene_count  = @total_genes;

            my $redisc_ratio      = ($redisc_gene_count / $orig_gene_count);
            my $gwide_ratio       = ($genomewide_count / $total_gene_count);

            my $enrichment        = ($redisc_ratio / $gwide_ratio);

            # format at last possible step:
            my $print_redisc_gene_count = commify($redisc_gene_count);
            my $print_orig_gene_count   = commify($orig_gene_count);

            my $print_genomewide_count  = commify($genomewide_count);
            my $print_total_gene_count  = commify($total_gene_count);

            my $print_redisc_ratio = sprintf "%.2f", $redisc_ratio;
            my $print_gwide_ratio  = sprintf "%.4f", $gwide_ratio;
            my $print_enrichment   = sprintf "%.2f", $enrichment;

            print "$header\n" if $header;
            $header = q{};

            print "$data_type\t";
            print "$motif_1\t";
            print "$print_redisc_gene_count\t";
            print "$print_orig_gene_count\t";
            print "$print_genomewide_count\t";
            print "$print_total_gene_count\t";
            print "$print_redisc_ratio\t";
            print "$print_gwide_ratio\t";
            print "$print_enrichment";
            print "\n";
        }
    }
    else { 
        die "Cannot parse input: $fimo_file\n";
    }    
}
close $LIST;


# Source -- Perl Cookbook 2.16, p. 84:
sub commify { 
    my $_text = reverse $_[0];
    $_text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $_text;
}

