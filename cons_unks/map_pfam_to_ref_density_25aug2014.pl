#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils qw(uniq);

my $data_ref;

my @pfams2uniprots = ();
my @uniprots2genes  = ();
my @genes2densities = ();

my $help;

GetOptions ( 'pfams2uniprots=s{,}'  => \@pfams2uniprots,
             'uniprots2genes=s{,}'  => \@uniprots2genes,
             'genes2densities=s{,}' => \@genes2densities,
             'help'                 => \$help, 
);


if ( $help or (! @pfams2uniprots) or (! @uniprots2genes) or (! @genes2densities) ) { 
    die "Format: map_pfam_to_ref_density.pl\n",
        "    --pfams2uniprots|-p   [1+ precomputed UniProt/PFAM tables from PFAM]\n",
        "    --uniprots2genes|-u   [1+ uniprot to gene tables; if using several, make sure no ambiguous short names!]\n",
        "    --genes2densities|-g  [1+ gene to annotation density tables; again, avoid short ambiguous names]\n",
        "    --help|-h             [print this message]\n",
        ;
}

foreach my $pfam2uniprot (@pfams2uniprots) { 
    open my $PFAM, '<', $pfam2uniprot or die "Can't open pfam2uniprot file $pfam2uniprot: $!";
    while (my $input = <$PFAM>) {
        chomp $input;
        # Sample input:
        # B5BM27	524	643	524	646	PF14843	GF_recep_IV	Domain	1	129	132	87.20	1e-21	CL0547
        if ( $input !~ /\A [#] /xms ) { 
            if ( $input =~ /\A (\S+) \t (?:[^\t]+ \t){4} (\S+) \t (\S+) \t /xms ) { 
                my $uniprot   = $1;
                my $pfam_id   = $2;
                my $pfam_text = $3;
                my $pfam_name = $pfam_id . q{|} . $pfam_text;
                $data_ref->{'pfam'}->{$pfam_name}->{'uniprot'}->{$uniprot} = 1;
            }
            else { 
                die "From pfam2uniprot file $pfam2uniprot, can't parse: $input\n";
            }
        }
    }
    close $PFAM or die "Can't close filehandle to pfam2uniprot file $pfam2uniprot: $!";
}

foreach my $uniprot2gene (@uniprots2genes) {
    open my $GENE, '<', $uniprot2gene or die "Can't open uniprot2gene file $uniprot2gene: $!";
    while (my $input = <$GENE>) {
        chomp $input;
        if ( $input =~ /\A (\S+) \t (\S+) \b/xms ) {
            my $uniprot = $1;
            my $gene    = $2;
            $data_ref->{'uniprot'}->{$uniprot}->{'gene'}->{$gene} = 1;
        }
    }
    close $GENE or die "Can't close filehandle to uniprot2gene file $uniprot2gene: $!";
}

foreach my $gene2density (@genes2densities) { 
    open my $DENSITY, '<', $gene2density or die "Can't open gene2density file $gene2density: $!";
    while (my $input = <$DENSITY>) {
        chomp $input;
        if ( $input =~ /\A (\S+) \t (\S+) \b/xms ) {
            my $gene    = $1;
            my $density = $2;
            $data_ref->{'gene'}->{$gene}->{'density'} = $density;
        }
    }
    close $DENSITY or die "Can't close filehandle to gene2density file $gene2density: $!";
}

my @pfam_domains = sort keys %{ $data_ref->{'pfam'} };
my @pfam_reports = ();
foreach my $pfam (@pfam_domains) { 
    my $total_ref_density = 0;
    my @genes             = ();
    my @uniprots = sort keys %{ $data_ref->{'pfam'}->{$pfam}->{'uniprot'} };
    foreach my $uniprot (@uniprots) {
        my @new_genes = sort keys %{ $data_ref->{'uniprot'}->{$uniprot}->{'gene'} };
        push @genes, @new_genes;
    }

    @genes         = sort @genes;
    @genes         = uniq @genes;
    my $gene_list  = join '; ', @genes;
    my $gene_no    = @genes;
    my $gene_count = q{};
    if ( $gene_no >= 2 ) { 
        $gene_count = "$gene_no genes";
    }
    if ( $gene_no == 1 ) {
        $gene_count = "$gene_no gene";
    }

    foreach my $gene (@genes) { 
        my $gene_ref_weight = 0;
        if ( exists $data_ref->{'gene'}->{$gene}->{'density'} ) { 
            $gene_ref_weight = $data_ref->{'gene'}->{$gene}->{'density'};
        }
        $total_ref_density = $total_ref_density + $gene_ref_weight;
        # No more than four decimal places, please!
        $total_ref_density = sprintf "%.4f", $total_ref_density;
    }
    push @pfam_reports, "$pfam\t$total_ref_density\t$gene_count\t$gene_list";
}

@pfam_reports = sort { rank_by_annot($b) <=> rank_by_annot($a) } @pfam_reports;

foreach my $pfam_report (@pfam_reports) { 
    print "$pfam_report\n";
}

sub rank_by_annot {
    my $_report = $_[0];
    if ( $_report =~ /\A \S+ \t (\S+) \b/xms ) { 
        my $_annot_density = $1;
        return $_annot_density;
    }
    else { 
        die "Subroutine rank_by_annot can't parse annotation report: $_report\n";
    }
}
