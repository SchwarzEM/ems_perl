#!/usr/bin/env perl

# prot2gene_table.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/2/2010.
# Purpose: given a set of proteome files, tabulate CDSes-to-genes, with none in >=2 proteomes; prerequisite for prot2gene_omcl.pl.
# 
# Note: this maps CDSes to unique genes; if it mapped *proteins*, the mapping could not be unique
#     because it is entirely possible for a single polypeptide to be encoded by >2 genes.

use strict;
use warnings;
use Getopt::Long;

# Import proteome files by '--argument [file list]'.
#     Each argument groups proteomes; each group is then assigned 
#     a parsing rule.  This allows modular addition of new parsing 
#     rules, as needed ad-hoc, within a stable frame of general, 
#     unchanging code.

my @wormgene = ();
my @augustus = ();
my @gene     = ();
my @simple   = ();
my @two_col  = ();
my @other    = ();

my $parse_wormgene = '\A > (\S+) \s .+ (WBGene\d+)';        # Official WB headers.
my $parse_augustus = '\A > ( (\S+) \. [t]{0,1} \d+)';       # AUGUSTUS-like headers (">gene.t1", ">gene.t2" or ">gene.1", ">gene.2").
my $parse_gene     = '\A > (\S+) \s .* \b [gG]ene: (\S+)';  # 'gene:' headers  (">CDS [...] [gG]ene:[gene_name]").
my $parse_simple   = '\A > (\S+) \s+ (\S+)';                # Simple headers,  (">CDS  gene").
my $parse_two_col  = '\A > (\S+) \s+ \S+ \s+ (\S+)';        # A new WormBase vogue in 3/2012 (">PROTEIN  [TRANSCRIPT?]  GENE").
my $parse_other    = '\A > ((\S+))';                        # Fall-back: each CDS name *is* the gene name.

my $data_ref;
my %cds2gene = ();

my $help;

GetOptions ( 'wormgene:s{,}' => \@wormgene,
             'augustus:s{,}' => \@augustus,
             'gene:s{,}'     => \@gene,
             'simple:s{,}'   => \@simple,
             'two_col:s{,}'  => \@two_col,
             'other:s{,}'    => \@other,
             'help'          => \$help,     );

my @infiles = (@wormgene, @augustus, @gene, @simple, @two_col, @other);

if ( ($help) or (! @infiles) ) { 
    die "Format: prot2gene_table.pl\n",
        "            --wormgene|-w (WBGene names)\n",
        "            --augustus|-a (AUGUSTUS-like, \"gene.t1+\" or \"gene.1+\")\n",
        "            --gene|-g     (\">CDS [...] [gG]ene:[gene_name]\")\n",
        "            --simple|-s   (\">CDS  [gene_name]\")\n",
        "            --two_col|-t  (\">CDS [cruft] [gene_name]\")\n",
        "            --other|-o    (\">CDS\" == gene name)\n",
        ;
}

foreach my $proteome (@wormgene) { 
    $data_ref->{'proteome'}->{$proteome}->{'rule'} = $parse_wormgene;
}
foreach my $proteome (@augustus) { 
    $data_ref->{'proteome'}->{$proteome}->{'rule'} = $parse_augustus;
}
foreach my $proteome (@gene) {
    $data_ref->{'proteome'}->{$proteome}->{'rule'} = $parse_gene;
}
foreach my $proteome (@simple) {
    $data_ref->{'proteome'}->{$proteome}->{'rule'} = $parse_simple;
}
foreach my $proteome (@two_col) { 
    $data_ref->{'proteome'}->{$proteome}->{'rule'} = $parse_two_col;
}
foreach my $proteome (@other) { 
    $data_ref->{'proteome'}->{$proteome}->{'rule'} = $parse_other;
}

foreach my $input_file (@infiles) { 
    my $rule = $data_ref->{'proteome'}->{$input_file}->{'rule'};
    open my $INFILE, '<', $input_file or die "Can't open input proteome file $input_file!\n";
    my %seen = ();
    while (my $input = <$INFILE>) { 
        chomp $input;

        # Only attempt to parse header lines.
        #    (If using non-FASTA-headers as input, must generalize this...)
        if ( $input =~ /\A >/xms ) { 

            # Enforce that the parsing rule always works (crucial!):
            if ( $input !~ /$rule/xms ) {
                warn "Input file $input_file is supposed to follow format:\n";
                warn "    ", $rule, "\n";
                die "Aberrant input: $input\n";
            }

            # Use validated rule to parse the FASTA header data:
            if ( $input =~ /$rule/xms ) { 
                my $cds  = $1;
                my $gene = $2;

                # Enforce only one CDS per input file:
                if (exists $seen{$cds} ) { 
                    die "CDS $cds occurs twice in input file $input_file!\n";
                }
                $seen{$cds} = 1;

                # Enforce unique CDS-to-gene throughout total of all input files:
                if ( ( exists $cds2gene{$cds} ) and ( $cds ne $cds2gene{$cds} ) ) {
                    die "CDS $cds is linked to both $cds2gene{$cds} and $gene!\n";
                }
                $cds2gene{$cds} = $gene;
            }
        }
    }
    close $INFILE or die "Can't close filehandle to proteome file $input_file!\n";
}

foreach my $cds (sort keys %cds2gene) { 
    print "$cds\t$cds2gene{$cds}\n";
}

