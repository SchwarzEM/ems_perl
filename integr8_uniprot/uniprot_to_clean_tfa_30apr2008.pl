#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;

# Program: uniprot_to_clean_tfa.pl -- Erich Schwarz <emsch@its.caltech.edu>, 4/30/2008.  LEGACY version, kept to ensure ability to reproduce past work.
# Purpose: Convert integr8.dat to FASTA with good headers.

# Have to switch between reading text and reading residues:
my $aa_read_toggle  = 0;
my $read_genes      = 0;
my $read_definition = 0;

# Initialize header variables:
my $protein     = q{}; 
my $uniprot_acc = q{};
my $gene        = 'no_gene_name';
my $definition  = q{};
my $species     = 'no_species_name';
my $ncbitax     = q{};
my $ncbi_acc    = q{};

while (my $input = <>) { 
    chomp $input;
    my $sourcefile = basename($ARGV);

    if ( $input =~ /\A ID \s+ (\S+) /xms ) {
        # 'ID' gives the canonical protein name.
        $protein        = $1;

        # DON'T try to get protein sequence yet.
        $aa_read_toggle = 0;

        # Start trying to get a gene name.
        $read_genes     = 1;

        # Initialize these to zero or boilerplate.
        $uniprot_acc    = q{};
        $gene           = 'no_gene_name';
        $species        = 'no_species_name';
        $ncbitax        = q{};

        # Later, script ignores all but the last-read NCBI acc.
        $ncbi_acc       = $protein;
    }

    elsif ( $input =~ /\A AC \s+ (\S+) ; /xms ) {
        $uniprot_acc = $1;
    }

    # Start to read multiline definitions to a single text line.
    elsif ( ( $input =~ /\A DE \s+ (\S.*) \s* \z/xms ) 
            && ( !$read_definition ) ) { 
        # Start recording a definition line:
        $definition = $1;
        $read_definition = 1;
    }
    # Stop reading when def. has a final  '.'!
    elsif ( ( $input =~ /\A DE \s+ (\S.*) \. \s* \z/xms )
            and ( $read_definition ) ) {
        $definition .= q{ } . $1;
        # Yes, stop:
        $read_definition = 0;
        # Condense catenated lines.
        $definition =~ s/[ ]{2,}/ /g;

        # Finally, store this line as 'gene':
        $gene = $definition;
        # Failsafe only -- keep looking for real gene name.
    }
    # If no period yet:
    elsif ( ( $input =~ /\A DE \s+ (\S.*) \s* \z/xms )
            and ( $read_definition ) ) {
        $definition .= q{ } . $1;
    }

    # The wonkiness of parsing the 'GN' line, pt. 1:
    elsif ( ( $input =~ / \A 
                          GN \s+ 
                          Names?= 
                          ([^;]+)    # $1
                          ; 
                          .+ 
                          Synonyms? 
                          = 
                          ([^;]+)    # $2 
                          ; 
                          .+ 
                          (?: OrderedLocusNames?= | ORFNames?=)
                          ([^;]+)    # $3
                          ; 
                       /xms )
            and ( $read_genes == 1 ) 
          ) { 
        my ($name1, $name2, $name3);
        ($name1, $name2, $name3) = ($1, $2, $3);

        # Consolidate names, incl. comma-ed list, to /-ed list.
        $gene = $name1 . '/' . $name2 . '/' . $name3;
        $gene =~ s/, /\//g;

        # Reformat gene name to one main, others parenthesesed.
        $gene =~ s#\A([^/]+)/(.+)\z#$1 \($2\)#;    # ?!

        # Stop trying to get a gene name.
        $read_genes = 0;
    }

    # The wonkiness of parsing the 'GN' line, pt. 2:
    elsif ( ( $input =~ / \A
                          GN\s+
                          Names?=
                          ([^;]+)  # $1
                          ; 
                          .+
                          (?: Synonyms?=|OrderedLocusNames?= | ORFNames?=)
                          ([^;]+)  # $2
                          ; 
                        /xms )
           and ( $read_genes == 1 )
          ) { 
        my ($name1, $name2);
        ($name1, $name2) = ($1, $2);

        # Consolidate names, incl. comma-ed list, to /-ed list.
        $gene = $name1 . '/' . $name2;
        $gene =~ s/, /\//g;

        # Reformat gene name to one main, others parenthesesed.
        $gene =~ s#\A([^/]+)/(.+)\z#$1 \($2\)#;

        # Stop trying to get a gene name.
        $read_genes = 0;
    }

    # The wonkiness of parsing the 'GN' line, pt. 3!:
    elsif ( (    ( $input =~ /\A GN.+ Names?=             ([^;]+) ;    /xms )
              or ( $input =~ /\A GN.+ Synonyms?=          ([^;]+) [,;] /xms )
              or ( $input =~ /\A GN.+ OrderedLocusNames?= ([^;]+) [,;] /xms )
              or ( $input =~ /\A GN.+ ORFNames?=          ([^;]+) [,;] /xms )
            )
        and ( $read_genes == 1 ) ) {
        $gene = $1;
        $gene =~ s/, /\//g;
        $gene =~ s#^([^/]+)/(.+)$#$1 \($2\)#;    # To cope with internal slashes.
        $read_genes = 0;
    }

    # Hopeless 'GN' line.  Annotate accordingly.
    elsif ( ( $input =~ /\A GN \s* /xms ) 
            and ( $read_genes == 1 ) ) {
        $gene       = 'confusing_gene_name';
        $read_genes = 0;
    }

    # And now! the 'DR' line:
    elsif ( ( $input =~ / \A
                          DR\s+
                          Ensembl;
                          \s+
                          (ENS[A-Z]+\d+) 
                          ; 
                          \s+
                          Homo sapiens 
                          .
                        /xms )
        and ( $read_genes == 1 ) ) {
        $gene       = $1;
        $read_genes = 0;
    }

    # Parsing the species line, pt. 1:
    elsif ( $input =~ / \A 
                        OS\s+ 
                        (.+) 
                        \. 
                        \s* 
                        \z 
                      /xms ) { 

        # What?
        if ( $species eq 'no_species_name' ) { 
            $species = q{}; 
        }
        # I have no idea why I coded this way, and it may be wrong:
        unless ( $species eq 'no_species_name' ) {
            $species .= $1;
        }
    }

    # Parsing the species line, pt. 2:
    elsif ( $input =~ / \A 
                        OS\s+ 
                        (.+[^.]) 
                        \s* 
                        \z 
                      /xms ) { 
        my $addendum = $1;

        # What?
        if ( $species eq "no_species_name" ) { 
            $species = q{}; 
        }
        # I have no idea why I coded this way, and it may be wrong:
        unless ( $species eq "no_species_name" ) { 
            $species = ( $species . $addendum . q{ } );
        }
    }

    # Extracting the NCBI taxon number:
    elsif ( $input =~ /\A OX\s+ (NCBI_TaxID=\d+) /xms ) {
        $ncbitax = $1;
    }

    # For some weird reason, it's the 'DR' line again:
    elsif ( $input =~ / DR\s+ 
                        EMBL 
                        ; 
                        \s+ 
                        \w+ 
                        ; 
                        \s+ 
                        ([A-Z]{3,}\d+\.\d+) 
                        ; 
                        .*? 
                        ; 
                        Genomic_DNA 
                        . 
                      /xms ) {
        $ncbi_acc = $1; 
        # N.B.: if multiple NCBI accs, repeatedly overwrites $ncbi_acc until final value.
    }

    # '//' == end of record.  Stop reading protein sequence.
    elsif ( ( $input =~ / \A \/ \/ /xms ) 
            and ( $aa_read_toggle == 1 ) ) {
        $aa_read_toggle = 0;
    }

    # Get ready to read protein seq., and print header pre-emptively.
    elsif ( $input =~ /\A SQ \s+ /xms ) {
        # Activate protein-reading:
        $aa_read_toggle = 1;

        # Printing a header line ...
        print ">$protein   $uniprot_acc   $ncbi_acc   $gene   $species ";

        # ... with bracketed NCBI tax. no. ...
        print q{[};
        print "$ncbitax";
        print q{]};

        # ... and the sourcefile named (where it's not just '-').
        print "   $sourcefile\n";
    }

    # Wow!  Actually reading protein sequence!
    elsif ( ( $input =~ /[a-zA-Z]+/xms ) 
            and ( $aa_read_toggle == 1 ) ) {
        my $orf_seq_line = $input;
        chomp($orf_seq_line);
        $orf_seq_line =~ s/[\d]*//g;
        $orf_seq_line =~ s/[\s]*//g;
        $orf_seq_line =~ s/[\W]*//g;
        $orf_seq_line =~ tr/a-z/A-Z/;
        print( $orf_seq_line . "\n" );
    }
}

