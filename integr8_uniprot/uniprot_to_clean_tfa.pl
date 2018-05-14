#!/usr/bin/env perl

use strict;
use warnings;
use autodie;
use Getopt::Long;
use File::Basename;

# Program: uniprot_to_clean_tfa.pl -- Erich Schwarz <emsch@its.caltech.edu>, 7/23/2012.
# Purpose: Convert EBI.dat (SwissProt or TREMBL) to FASTA with good headers, optionally, prepend species_specific prefixes to sequence names.

my @input_files = ();

my $prefix = q{};
my $help;

# To be used as filehandle later:
my $INPUT_FILE;

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

GetOptions ( 'input_files=s{,}' => \@input_files,
             'prefix=s'         => \$prefix,
             'help'             => \$help, );

if ( $help or (! @input_files) ) { 
    die "\n",
        "Format: uniprot_to_clean_tfa.pl\n",
        "        --input_files|-i  [input file(s), or '-' if stream]\n",
        "        --prefix|-p       [make species-specific prefixes, with default prefix named here; prefix will always have a single final '_' character]\n", 
        "        --help|-h         [print this message]\n",
        "\n",
        "Example: uniprot_to_clean_tfa.pl -i swissprot.dat trembl.dat -p EMBL [prepended as 'EMBL_', *if* no species prefix successfully made] > monster.fa ;\n",
        "\n",
        ;
}

# *If* we have decided to prepend a prefix, always give it a single final '_' character, 
#     silently censoring any *additional* '_' characters at the end.
# Otherwise, leave it blank!
if ($prefix) { 
    $prefix =~ s/[_]+\z//;
    $prefix = $prefix . q{_};
}

# Accept either a stream from '-' or a standard file.
foreach my $infile (@input_files) { 
    if ($infile eq '-') {
        # Special case: get the stdin handle
        $INPUT_FILE = *STDIN{IO};
    }
    else {
        # Standard case: open the file
        open $INPUT_FILE, '<', $infile;
    }

    # Record the incoming FASTA data.
    while (my $input = <$INPUT_FILE>) {

        chomp $input;
        my $sourcefile = basename($infile);

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

            # Printing a header line.

            # This used to space the fields with '   '.
            # But it now uses "\t" between fields, so that the header can be easily subjected to further parsing, as needed.
            # Note that we are now opting to treat $species and $ncbitax as different fields, and that we are not putting $ncbitax in [...];
            #     these are both significant changes from the previous 2008 version of uniprot_to_clean_tfa.pl.
            # The sourcefile is named, where it's not just '-'.

            # If we have not asked for a prefix, $prefix will be q{}, and thus be invisible here.
            # If we have asked for a prefix, we will get a species-specific prefix if possible; otherwise we will get a default.

            # We keep $prefix as a stable default; but we individually specify prefixes for each individual sequence.
            my $indiv_prefix = $prefix;
            if ($prefix) { 
                if ( $species =~ /\A (\S+) \s+ (\S+) /xms ) { 
                    my $first_name  = $1;
                    my $second_name = $2;
                    my $len_first_name  = length $first_name;
                    my $len_second_name = length $second_name;
                    if ( ( $len_first_name >= 2 ) and ( $len_second_name >= 3) ) {
                        $first_name  =~ s/\A(\S{2})\S*\z/$1/;
                        $second_name =~ s/\A(\S{3})\S*\z/$1/;
                    }
                    elsif ( ( $len_first_name >= 3 ) and ( $len_second_name >= 2) ) {
                        $first_name  =~ s/\A(\S{3})\S*\z/$1/;
                        $second_name =~ s/\A(\S{2})\S*\z/$1/;
                    }
                    elsif ( ( $len_first_name >= 4 ) and ( $len_second_name >= 1) ) {
                        $first_name  =~ s/\A(\S{4})\S*\z/$1/;
                        $second_name =~ s/\A(\S{1})\S*\z/$1/;
                    }
                    $indiv_prefix = $first_name . $second_name . q{_};
                }
                elsif ( $species =~ /\A (\S+) /xms ) {
                    my $first_name  = $1;
                    my $len_first_name  = length $first_name;
                    if ( $len_first_name >= 6 ) { 
                        $first_name  =~ s/\A(\S{5})\S*\z/$1/;
                    }
                    if ( $len_first_name >= 3 ) {
                        $indiv_prefix = $first_name . q{_};
                    }
                    # If we can't even get 3 letter from the first name, we fall back on the default prefix.
                }
                # If both of these fail, then we should get the default, user-specified prefix.
            }
            print ">$indiv_prefix$protein\t$uniprot_acc\t$ncbi_acc\t$gene\t$species\t$ncbitax\t$sourcefile\n";
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
    close $INPUT_FILE;
}

