#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils qw(uniq);
use autodie;

my @acey_data = ();
my @cel_data  = ();

my $data_ref;
my $help;

GetOptions ( 'acey=s{,}' => \@acey_data,
             'cel=s{,}'  => \@cel_data,
             'help'      => \$help,     );

if ( $help or (! @acey_data ) or (! @cel_data ) ) { 
    die "Format: make_pre_drug_summary_21oct2014.pl\n",
        "    --acey|-a   [N files, with drug target data for A. ceylanicum genes]\n",
        "    --cel|-c    [N files, with drug target data for C. elegans genes]\n",
        "    --help|-h   [print this message]\n",
        ;
}

foreach my $acey_file (@acey_data) { 
    &read_file($acey_file, 'acey');
}

foreach my $cel_file (@cel_data) { 
    &read_file($cel_file, 'cel');
}

my @acey_genes = sort keys %{ $data_ref->{'species'}->{'acey'}->{'gene'} };
foreach my $acey_gene (@acey_genes) {
    my @acey_motifs          = sort keys %{ $data_ref->{'species'}->{'acey'}->{'gene'}->{$acey_gene}->{'motif'} };
    my @acey_motif_annots    = ();
    my @acey_motif_cel_genes = ();
    foreach my $acey_motif (@acey_motifs) {
        my @acey_motif_genes   = sort keys %{ $data_ref->{'species'}->{'acey'}->{'motif'}->{$acey_motif}->{'gene'} };        
        my $acey_motif_gene_no = @acey_motif_genes;
        my $acey_motif_annot   = "$acey_motif ($acey_motif_gene_no genes)";
        push @acey_motif_annots, $acey_motif_annot;
        my @mot_assoc_cel_genes = sort keys %{ $data_ref->{'species'}->{'cel'}->{'motif'}->{$acey_motif}->{'gene'} };
        push @acey_motif_cel_genes, @mot_assoc_cel_genes;
    }
    @acey_motif_annots    = sort @acey_motif_annots;
    @acey_motif_cel_genes = sort @acey_motif_cel_genes;
    @acey_motif_annots    = uniq @acey_motif_annots;
    @acey_motif_cel_genes = uniq @acey_motif_cel_genes;
    my $acey_motif_txt    = join '; ', @acey_motif_annots;
    my $acey_cel_gene_txt = join '; ', @acey_motif_cel_genes;
    print "$acey_gene\t$acey_motif_txt\t$acey_cel_gene_txt\n";
}

sub read_file {
    my $input_file = $_[0];
    my $species    = $_[1];
    my $not_header = 0;
    open my $FILE, '<', $input_file;
    while (my $input = <$FILE>) {
        # Parse every line of the table, *except* the header line.
        if ($not_header) {          
            chomp $input;
            my $gene        = q{};
            my $motif_title = q{};
            # For PFAM table:
            if ( $input =~ /\A ([^\t]+) \t ([^\t]+) \t ([^\t]+) \z/xms ) {
                my $motif_acc   = $1;
                my $motif_name  = $2;
                $gene           = $3;
                $motif_title    = "$motif_name [$motif_acc]";
            }
            # For InterPro table:
            elsif ( $input =~ /\A ([^\t]+) \t ([^\t]+) \t ([^\t]+) \t ([^\t]+) \z/xms ) {
                my $motif_type  = $1;
                my $motif_acc   = $2;
                my $motif_name  = $3;
                $gene           = $4;
                $motif_title    = "$motif_name [" . $motif_type . q{|}. "$motif_acc]";
            }
            else { 
                die "From input file $input_file, can't parse input line: $input\n";
            }
            $data_ref->{'species'}->{$species}->{'gene'}->{$gene}->{'motif'}->{$motif_title} = 1;
            $data_ref->{'species'}->{$species}->{'motif'}->{$motif_title}->{'gene'}->{$gene} = 1;
        }
        else {
            $not_header = 1;
        }
    }
    close $FILE;
    return;
}

