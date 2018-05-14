#!/usr/bin/env perl

# add_tab_annots_unreliable_legacy.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/22/2010.
# Purpose: given index genefile and 1+ nonredund. annot. TSVs, append annots to index text by "\t" or "; ".  Note that this is a defective script which fails to properly fill in blank tabs for genes lacking annotations in a particular file!  It is being kept around for purely legacy reasons (i.e., to be able to reproduce my own awful code if need be).  Use add_tab_annots.pl instead.

use strict;
use warnings;
use Getopt::Long;

my $index = q{};
my @annot_files = ();
my $full_annot;
my $semicolon;
my $help;

my $gene         = q{};
my $annot        = q{};
my $i            = 0;
my $single_field = q{};

my $index_lines_ref;

GetOptions ( 'index=s'     => \$index,
             'annots=s{,}' => \@annot_files,
             'full_annot'  => \$full_annot,
             'semicolon'   => \$semicolon,
             'help'        => \$help, );

if ($help or (! $index) or (! @annot_files) ) { 
    die "Format: add_tab_annots_unreliable_legacy.pl [why are you still using this?]\n",
        "    --index|-i [index file to have annotations appended to it]\n",
        "    --full_annot|-f [use full text on annotation line, not just last field]\n",
        "    --semicolon|-s [replace ", q{"." with "; "}, " default is to append ", q{"\t"}, "]\n",
        "    --annots|-a <list of annot files; nonredundant annots./file>\n",
        ;
}

open my $INDEX, '<', $index or die "Can't open index file $index: $!";
while (my $input = <$INDEX>) { 
    chomp $input;
    if ($input =~ /\A (\S+) \b /xms ) { 
        $gene = $1;
        $index_lines_ref->{'i'}->{$i} = $gene;
        if ( exists $index_lines_ref->{'gene'}->{$gene} ) { 
            die "Gene $gene in index file $index has redundant annotation: $input\n";
        }
        $index_lines_ref->{'gene'}->{$gene}->{'text'} = $input;
        $i++;
    }
    else { 
        die "Can't parse input from index file $index: $input\n";
    }
}
close $INDEX or die "Can't close filehandle to index file $index: $!";

foreach my $annot_file (@annot_files) { 
    if ( $annot_file !~ /\w/ ) {
        die "Unusable annotation file name: $annot_file\n";
    }
    open my $AFILE, '<', $annot_file or die "Can't open annotation file $annot_file: $!";
    while (my $input = <$AFILE>) { 
        chomp $input;
        # The following was '(.+)', but I decided that this was bad, because it kept me from adding empty columns.
        if ($input =~ /\A ([^\t]+) \t (.*) \z /xms ) { 
            $gene  = $1;
            $annot = $2;
            # The *default* is to only read the last tab-delimited data field:
            if ( (! $full_annot ) and ( $annot =~ / \A .* \t ([^\t]+) \z /xms ) ) { 
                $single_field = $1;
                $annot = $single_field;
            }
            if ( exists $index_lines_ref->{'gene'}->{$gene}->{'text'} ) { 
                if ( exists $index_lines_ref->{'gene'}->{$gene}->{'annot'}->{$annot_file} ) { 
                    die "Redundant gene annotation from $annot_file in: $input\n";
                }
                $index_lines_ref->{'gene'}->{$gene}->{'annot'}->{$annot_file} = $annot;
            }
        }
        else { 
            die "Can't parse input from added annotation file $annot_file: $input\n";
        }
    }
    close $AFILE or die "Can't close filehandle to annotation file $annot_file: $!";
}

# Append annotations in the order specified for input files (if any).
# Remember that we counted lines from 0 to $i-1, and never did actually append a "line $i" for the last line.

my $j = $i;
$j--;

if ( exists $index_lines_ref->{'i'}->{$i} ) {
    die "Should not have a data field for i = $i\n";
}

foreach my $line_no (0 .. $j) { 
    if ( exists $index_lines_ref->{'i'}->{$line_no} ) { 
        my $annot_gene = $index_lines_ref->{'i'}->{$line_no};
        if (! $semicolon) { 
            print "$index_lines_ref->{'gene'}->{$annot_gene}->{'text'}";
            foreach my $append_file (@annot_files) { 
                print "\t";
                if (exists $index_lines_ref->{'gene'}->{$annot_gene}->{'annot'}->{$append_file} ) { 
                    print "$index_lines_ref->{'gene'}->{$annot_gene}->{'annot'}->{$append_file}";
                }
            }
            print "\n";
        }
        if ($semicolon) { 
            if ($full_annot) { 
                die "Not yet able to handle both semicolons and complex annotations!\n";
            }
            my $text = $index_lines_ref->{'gene'}->{$annot_gene}->{'text'};
            my $endtext = q{};
            if ( $text =~ / \A (.*?) ([.]{0,1} ["]{0,1}) \s* \z /xms ) { 
                $text    = $1;
                $endtext = $2;
                print "$text";
                foreach my $append_file (@annot_files) { 
                    if (exists $index_lines_ref->{'gene'}->{$annot_gene}->{'annot'}->{$append_file} ) { 
                        print "; $index_lines_ref->{'gene'}->{$annot_gene}->{'annot'}->{$append_file}";
                    }
                }
                print "$endtext\n";
            }
        }
    }
    else { 
        warn "Failed to parse concatenated data for line $j:\n";
    }
}

