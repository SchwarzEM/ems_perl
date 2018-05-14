#!/usr/bin/env perl

# genes_from_taxon_in_omcl.pl -- Erich Schwarz <emsch@its.caltech.edu>, 10/22/2009; upgraded 11/2/2012; tweaked to allow OrthoFinder input in OrthoMCL-ized format, 1/31/2017.
# Purpose: given OrthoMCL report text (file or piped) and a named taxon, count genes/proteins of that taxon and summarize them; optionally, provide only the summary or the list.

use strict;
use warnings;
use Getopt::Long;

my @infiles       = ();
my $taxon         = q{};
my $genes_only;
my $summary_only;
my $help;

my $group_count   = 0;
my $pos_grp_count = 0;
my $gene_count    = 0;

my %genes = ();

GetOptions ( 'infile=s{,}' => \@infiles,
             'taxon=s'     => \$taxon,
             'genes'       => \$genes_only,
             'summary'     => \$summary_only,
             'help'        => \$help,  );

if ( $help or (! $taxon ) or ( $genes_only and $summary_only ) ) {
    die "Format: genes_from_taxon_in_omcl.pl\n",
        "        --taxon|-t    [taxon to count]\n", 
        "        --infile|-i   [input files or '-' for stream]\n",
        "        --genes|-g    [only list taxon-specific genes/proteins in the OrthoMCL report; mutually exclusive with --summary]\n",
        "        --summary|-s  [only summarize the count of taxon-specific genes/proteins; default is to print both summary and gene list]\n",
        "        --help|-h     [print this message]\n",
        ;
}

# To be used as filehandle for the input(s):
my $INFILE;

foreach my $infile (@infiles) { 
    if ($infile eq '-') {
        # Special case: get the stdin handle
        $INFILE = *STDIN{IO};
    }
    else {
        # Standard case: open the file
        open $INFILE, '<', $infile or die "Can't open input file $infile. $!\n";
    }
    while (my $input = <$INFILE>) { 
        chomp $input;

        # Sample input line:  
        # ORTHOMCL19007(2 genes,2 taxa: WBGene00031156(briggsae) WBGene00005532(elegans)

        if ( ( $input =~ /\S/xms ) 
             and ( $input !~ / \A 
                               [A-Z]+\d+    # revise this to allow OrthoFinder results (reworked into OrthoMCL format) 
                               \( \d+ \s+ genes,\d+ \s+ taxa\) : 
                               (?: \s+ [^\s\(\)]+ \( [^\s\(\)]+ \) )+
                               \s* \z /xms ) ) { 
            die "Obviously not an OrthoMCL text line!\n";
        }
        else { 
            $group_count++;
            if ( $input =~ / \( $taxon \) /xms ) { 
                $pos_grp_count++;
            }
            while ( $input =~ / (\S+) \( $taxon \) /gxms ) { 
                my $gene = $1;
                $genes{$gene} = 1;
                $gene_count++;
            }
        }
    }
    close $INFILE or die "Can't close filehandle to input file $infile. $!\n";
}

if (! $genes_only ) { 
    $group_count   = commify($group_count);
    $pos_grp_count = commify($pos_grp_count);
    $gene_count    = commify($gene_count);

    print "\n";
    print "                Taxon: $taxon\n";
    print "         Total groups: $group_count\n";
    print "Taxon-positive groups: $pos_grp_count\n";
    print "     Genes from taxon: $gene_count\n";
    print "\n";
}

if ( (! $genes_only ) and (! $summary_only ) ) { 
    print "List of taxon-specific genes:\n";
    print "----------------------------\n";
}

if (! $summary_only ) { 
    my @gene_list = sort keys %genes;
    foreach my $gene ( @gene_list ) { 
        print "$gene\n";
    }
}

sub commify { 
    my $text = reverse $_[0];
    $text =~ s/ (\d{3}) 
                 (?=\d) 
                 (?!\d*\.)
               /$1,/xmsg;
    return scalar reverse $text;
}

