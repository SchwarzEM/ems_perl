#!/usr/bin/env perl

# i8fa_gp_2sp_pr_roster.pl -- Erich Schwarz <emsch@its.caltech.edu>, 5/12/2008.
# Purpose: given integr8-der. FASTA, UniProt->RefSeq tab., and GenPept, make species + (proteins) list.

use strict;
use warnings;
use Getopt::Long;

my $fasta             = q{};
my $refs_file         = q{};
my $gpep_file         = q{};
my $taxonline         = q{};
my %taxa2uniprots_ref = ();

GetOptions ( "fasta=s"  => \$fasta,  
             "refseq=s" => \$refs_file, 
             "genpep=s" => \$gpep_file, );

my %seen_in_fasta     = ();
my %refseq2uprot      = ();
my %gpep_prots        = ();
my %refseq2taxon      = ();
# my %species2rseqs_ref = ();

if (    (! -r $fasta) 
     or (! -r $refs_file) 
     or (! -r $gpep_file ) ) { 
    die 'Format:',
        ' ./i8fa_gp_2sp_pr_roster.pl',
        ' --fasta=[intgr8-d. FASTA]',
        ' --refseq=[UProt->RefSq tab.]',
        ' --genpep=[FASTA->RefSq-d. GenPept]',
        "\n",
       ;
}

open my $FASTA, '<', $fasta 
    or die "Can't open FASTA file $fasta: $!";
# Sample input:
# >A0AHB4_LISW6   A0AHB4 ... [etc.]
while (my $input = <$FASTA>) { 
    chomp $input;
    if ( $input =~ / \A > ((\w+)_\w+) /xms) { 
        my ($protname, $uniprot_id) = ($1, $2);
        $seen_in_fasta{$protname} = $uniprot_id;
        $seen_in_fasta{$uniprot_id} = $protname;
    }
}
close $FASTA or die "Can't close filehandle of $fasta: $!";

open my $UPROT_2_REFSEQ, '<', $refs_file 
    or die "Can't open UniProt->RefSeq table $refs_file: $!";
# Sample input:
# A0AHB4_LISW6	YP_849177.1
while (my $input = <$UPROT_2_REFSEQ>) {
    chomp $input;
    if ( $input =~ / \A (\S+_\S+) \s+ (\S+_\S+) /xms) {
        my ($uniprot_id, $refseq_id) = ($1, $2);
        if ($seen_in_fasta{$uniprot_id}) { 
            $refseq_id =~ s/\.\d+//;
            $refseq2uprot{$refseq_id} = $uniprot_id;
        }
    }
}
close $UPROT_2_REFSEQ 
    or die "Can't close filehandle of $refs_file: $!";

parse_genpept($gpep_file);

# Mostly made a sub to encapsulate passed-around $refseq_id...
sub parse_genpept { 
    my $refseq_id    = q{};
    my $species      = q{};
    my $taxonline    = q{};
    my $read_species = 'no';
    my $read_taxa    = 'no';


    open my $GENPEPT, '<', $gpep_file 
        or die "Can't open GenPept $gpep_file: $!";

# Sample start input line:
# LOCUS       YP_849177                343 aa            linear   BCT 27-MAR-2008

# Sample taxon input:
#   
#   ORGANISM  Mesorhizobium loti MAFF303099
#             Bacteria; Proteobacteria; Alphaproteobacteria; Rhizobiales;
#             Phyllobacteriaceae; Mesorhizobium.


    while (my $input = <$GENPEPT>) { 
        chomp $input;

        # An if-elsif-elsif chain is less than elegant code, but
        #    proved necessary to untangle parsing logic.

        if ( $input =~ / \A LOCUS \s+ (\S+_\S+) /xms) { 
            $refseq_id    = $1;
            $species      = q{};
            $taxonline    = q{};
            $read_species = 'no';
            $read_taxa    = 'no';
            if ( $refseq2uprot{$refseq_id} ) { 
                $gpep_prots{$refseq_id} = $refseq2uprot{$refseq_id};
            }
        }

        elsif ( ( $read_species eq 'no' ) 
              and ( $read_taxa eq 'no' ) 
              and ( $input =~ / \A \s+ ORGANISM \s+ (\S.*\S) \s* \z /xms ) ) { 
            $species = $1;
            $read_species = 'yes';
            $read_taxa    = 'no';
        }

        elsif ( ( $read_species eq 'yes' ) and ( $read_taxa eq 'no' ) ) { 
            # Absence of "; " means a confusing multiline species name...
            if ( ( $input !~ /;[ ]/xms) and ( $input =~ /\w/xms ) ) { 
                $input =~ s/\A\s*//;
                $input =~ s/\s*\z//;
                $species .= q{ } . $input;
            }

            # Taxon lines *do* have "; ", so switch if you see this.
            if ( $input =~ /;[ ]/xms) { 
                $read_species = 'no';
                $read_taxa = 'yes';
                $input =~ s/\A\s*//;
                $input =~ s/\s*\z//;
                # Start a new taxonline:
                $taxonline = $input;
            }
        }

        elsif ( ( $read_taxa eq 'yes' ) and ($read_species eq 'no' ) ) { 
            # If there's a new header, stop reading and finish taxonline.
            if ( $input =~ / \A [A-Z]+ /xms ) { 
                $read_species = 'no';
                $read_taxa = 'no';
                $taxonline =~ s/\.\z//;
                $taxonline .= '; ' . $species;
                $refseq2taxon{$refseq_id} = $taxonline;
                $species   = q{};
                $taxonline = q{};
            }
            else { 
                $input =~ s/\A\s*//;
                $input =~ s/\s*\z//;
                $taxonline .= q{ } . $input;
            }
        }
    }
    close $GENPEPT 
        or die "Can't close filehandle of $gpep_file: $!";
}
# End &parse_genpept().

my %output_lines = ();

foreach my $prot (sort keys %gpep_prots) { 
    my $taxon = $refseq2taxon{$prot};
    push @{ $taxa2uniprots_ref{$taxon} }, $refseq2uprot{$prot};
}

foreach my $taxon (sort keys %taxa2uniprots_ref) { 
    print $taxon,
          "\t",
          "[",
          "@{ $taxa2uniprots_ref{$taxon} }",
          "]",
          "\n",
          ;
}

