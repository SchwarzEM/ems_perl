#!/usr/bin/env perl

# phenoace2annot.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/22/2010.
# Purpose: given a .ace file of alleles or RNAis, extract gene-phenotype table or gene-NotPhenotype table.

use strict;
use warnings;
use Getopt::Long;

# Starting arguments:
my $phenonames = q{};
my $help;

# Data structure for mapping phenotype numbers to readable names:
my $wbpheno2text_ref;

# Variables for regex capture:
my $wbgene    = q{};
my $wbpheno   = q{};
my $phenotext = q{};

# Variables for tracking all phenotypes affecting all genes, in a single allele/RNAi.
my %genes           = ();
my %positive_phenos = ();   # Actual phenotypes seen.
my %negative_phenos = ();   # Phenotype scored but *not* seen ('Not').

# Data structure for storing final mappings of genes to positive/negative phenotypes.
my $gene_info_ref;

GetOptions ( 'phenonames=s' => \$phenonames,
             'help'         => \$help, );

if ($help or (! $phenonames) or (! @ARGV) ) { 
    die "Format: phenoace2annot.pl --phenonames|-p [.ace of phenotype names] <input .ace file(s) or stream>\n";
}

# Get data required to replace WBPhenotype:[number] with readable descriptions.
open my $PHENONAMES, '<', $phenonames or die "Can't open .ace file of phenotype names $phenonames: $!";
while (my $input = <$PHENONAMES>) { 
    chomp $input;
    if ( $input =~ / \A Phenotype \s+ : \s+ \" WBPhenotype: (\d+) \" /xms ) { 
        $wbpheno = $1;
        # This is to ensure that:
        #     (1) we always have a mapping, and
        #     (2) we catch failures to extract human-readable text.
        $wbpheno2text_ref->{$wbpheno}->{'seen'} = 1;
    }
    elsif ( $input =~ / \A Primary_name \s+ \" ([^\"]+) \" /xms ) { 
        $phenotext = $1;
        $wbpheno2text_ref->{$wbpheno}->{'primary'} = $phenotext;
    }
    elsif ( $input =~ / \A Short_name \s+ \" ([^\"]+) \" /xms ) {  
        $phenotext = $1;
        $wbpheno2text_ref->{$wbpheno}->{'short'} = $phenotext;
    }   
}
close $PHENONAMES or die "Can't close filehandle of .ace file of phenotype names $phenonames: $!";

foreach my $wbpheno2 ( sort keys %{ $wbpheno2text_ref } ) { 
    if ( exists $wbpheno2text_ref->{$wbpheno2}->{'short'} ) { 
        $wbpheno2text_ref->{$wbpheno2}->{'name'} 
            = $wbpheno2text_ref->{$wbpheno2}->{'short'};
    }
    elsif ( exists $wbpheno2text_ref->{$wbpheno2}->{'primary'} ) {
        $wbpheno2text_ref->{$wbpheno2}->{'name'}
            = $wbpheno2text_ref->{$wbpheno2}->{'primary'};
    }
    else { 
        $wbpheno2text_ref->{$wbpheno2}->{'name'} = $wbpheno2;
    }
}

while (my $input = <>) { 
    chomp $input;
    # Sample header lines:
    # RNAi : "WBRNAi00000001"
    # Variation : "a83"
    if ( $input =~ / \A (?:RNAi|Variation) \s+ : \s+ \"[^\"]+ \" /xms ) { 
        # Clear out stored data from previous variant:
        foreach my $g (sort keys %genes) { 
            foreach my $pphen (sort keys %positive_phenos) { 
                $gene_info_ref->{$g}->{'pos_phenos'}->{$pphen} = 1;
            }
            foreach my $nphen (sort keys %negative_phenos) { 
                $gene_info_ref->{$g}->{'neg_phenos'}->{$nphen} = 1;
            }
        }
        %genes           = ();
        %positive_phenos = ();
        %negative_phenos = ();
        $wbgene          = q{};
        $wbpheno         = q{};
    }
    elsif ( $input =~ / \A Gene \s+ \" (WBGene\d+) \" /xms ) { 
        $wbgene = $1;
        $genes{$wbgene} = 1;
        # To prevent inadvertant carry-over of scratch variables.
        $wbgene = q{};
    }
    elsif ( $input =~ / \A Phenotype \s+ \"WBPhenotype:(\d+)\" /xms ) { 
        $wbpheno = $1;
        # Map from serial numbers to human-readable texts:
        $wbpheno = $wbpheno2text_ref->{$wbpheno}->{'name'};
        if ( $input =~ / \A Phenotype \s+ \"WBPhenotype:(\d+)\" \s+ Not \b /xms ) { 
            $negative_phenos{$wbpheno} = 1;
        }
        else { 
            $positive_phenos{$wbpheno} = 1;
        }
        # To prevent inadvertant carry-over of scratch variables.
        $wbpheno = q{};
    }
    elsif ( $input =~ / \A Phenotype_not_observed \s+ \"WBPhenotype:(\d+)\" /xms ) {
        $wbpheno = $1;
        # Again, map from serial numbers to human-readable texts:   
        $wbpheno = $wbpheno2text_ref->{$wbpheno}->{'name'};
        $negative_phenos{$wbpheno} = 1;
        $wbpheno = q{};
    }
}

# End-of-loop: do last clearing out of stored data from previous variant:
foreach my $g (sort keys %genes) {
    foreach my $pphen (sort keys %positive_phenos) {
        $gene_info_ref->{$g}->{'pos_phenos'}->{$pphen} = 1;   
    }
    foreach my $nphen (sort keys %negative_phenos) {
        $gene_info_ref->{$g}->{'neg_phenos'}->{$nphen} = 1;
    }
}
# Pro-forma:
%genes           = ();
%positive_phenos = ();
%negative_phenos = ();
$wbgene          = q{};
$wbpheno         = q{};

foreach my $gid (sort keys %{ $gene_info_ref }) { 

    my @positive_data  = sort keys %{ $gene_info_ref->{$gid}->{'pos_phenos'} };
    my $positive_annot = join "; ", @positive_data;
    if ($positive_annot) {
        $positive_annot = q{"WBPhenotype: } . $positive_annot . q{"};
    }

    my @negative_data  = sort keys %{ $gene_info_ref->{$gid}->{'neg_phenos'} }; 
    @negative_data = grep { ! exists $gene_info_ref->{$gid}->{'pos_phenos'}->{$_} } @negative_data;
    my $negative_annot = join "; ", @negative_data;
    if ($negative_annot) { 
        $negative_annot = q{"NOT WBPheno: } . $negative_annot . q{"};
    }

    print "$gid\t$positive_annot\t$negative_annot\n";
}

