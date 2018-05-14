#!/usr/bin/env perl

# add_tab_annots.pl -- Erich Schwarz <ems394@cornell.edu>, 12/10/2017.
# Purpose: given index file (of genes or other key entities) and 1+ nonredund. annot. TSVs, append annots to index text by "\t" or "; ".  Note that this is better engineered for tabbed input than semicolon input; it will enforce uniform tabcounts for each input file, then reliably fills in tabs for genes (or other key IDs) that do not have annotations.  This needs to also be implemented for semicolon input, but is considerably less important, so gets short shrift.

use strict;
use warnings;
use Getopt::Long;

my $index       = q{};
my @annot_files = ();
my $key_id      = 'Gene';  # default value

my $full_annot;
my $semicolon;
my $help;

my $annot        = q{};
my $i            = 0;

my $tab_count    = 0;

my $index_lines_ref;

GetOptions ( 'index=s'     => \$index,
             'annots=s{,}' => \@annot_files,
             'key_id'      => \$key_id,
             'full_annot'  => \$full_annot,
             'semicolon'   => \$semicolon,
             'help'        => \$help, );

if ($help or (! $index) or (! @annot_files) ) { 
    die "Format: add_tab_annots.pl\n",
        "    --index|-i [index file to have annotations appended to it]\n",
        "    --full_annot|-f [use full text on annotation line, not just last field]\n",
        "    --key_id|-k [Key ID / column name of first data field; defaults to \"Gene\"]\n",
        "    --semicolon|-s [replace ", q{"." with "; "}, " default is to append ", q{"\t"}, "]\n",
        "    --annots|-a <list of annot files; nonredundant annots./file>\n",
        ;
}

open my $INDEX, '<', $index or die "Can't open index file $index: $!";
$tab_count = 0;
while (my $input = <$INDEX>) { 
    chomp $input;
    if ($input =~ /\A (\S+) \b /xms ) { 
        $key_id = $1;
        $index_lines_ref->{'i'}->{$i} = $key_id;
        if ( exists $index_lines_ref->{'key_id'}->{$key_id} ) { 
            die "Key ID (e.g., gene) $key_id in index file $index has redundant annotation: $input\n";
        }
        $index_lines_ref->{'key_id'}->{$key_id}->{'text'} = $input;
        $i++;

        if (! $semicolon) {
            # Forbid ragged tabbed input index data.
            my $line_tab_count = ( $input =~ tr/\t/\t/ );
            if ( ($tab_count) and ( $line_tab_count != $tab_count ) ) { 
                die "Index file $index has inconsistent tab count in input lines: $tab_count vs. $line_tab_count\n";
            }
            $tab_count = $line_tab_count;
        }
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
    $tab_count = 0;
    while (my $input = <$AFILE>) { 
        chomp $input;
        # The following was '(.+)', but I decided that this was bad, because it kept me from adding empty columns.
        if ($input =~ /\A ([^\t]+) \t (.*) \z /xms ) { 
            $key_id  = $1;
            $annot = $2;
            # The *default* is to only read the last tab-delimited data field:
            if ( (! $full_annot ) and ( $annot =~ / \A .* \t ([^\t]+) \z /xms ) ) { 
                my $single_field = $1;
                $annot = $single_field;
            }
            if ( exists $index_lines_ref->{'key_id'}->{$key_id}->{'text'} ) { 
                if ( exists $index_lines_ref->{'key_id'}->{$key_id}->{'annot'}->{$annot_file} ) { 
                    die "Redundant key ID (e.g., gene) annotation from $annot_file in: $input\n";
                }
                $index_lines_ref->{'key_id'}->{$key_id}->{'annot'}->{$annot_file} = $annot;

                if (! $semicolon) {
                    # Again, forbid ragged tabbed input data -- in $annot, not $input!
                    my $line_tab_count = ( $annot =~ tr/\t/\t/ );
                    if ( ($tab_count) and ( $line_tab_count != $tab_count ) ) {
                        die "Annotation file $annot_file has inconsistent tab count in input lines: $tab_count vs. $line_tab_count\n";
                    }
                    $tab_count = $line_tab_count;
                    $index_lines_ref->{'annot_file'}->{$annot_file}->{'tab_count'} = $tab_count;
                }
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
        my $annot_key_id = $index_lines_ref->{'i'}->{$line_no};
        if (! $semicolon) { 
            print "$index_lines_ref->{'key_id'}->{$annot_key_id}->{'text'}";
            foreach my $append_file (@annot_files) { 
                # This is always needed as a spacer between columns of appended data:
                print "\t";

                # If we actually have an annotation, print that:
                if (exists $index_lines_ref->{'key_id'}->{$annot_key_id}->{'annot'}->{$append_file} ) { 
                    print "$index_lines_ref->{'key_id'}->{$annot_key_id}->{'annot'}->{$append_file}";
                }

                # Conversely, if we are using tabbed data, and we *don't* have an annotation:
                #     generate a blank line of tabs, which can be used safely to fill in empty annotations;
                #     then, print the empty tabs instead.

                if ( (! $semicolon) and (! exists $index_lines_ref->{'key_id'}->{$annot_key_id}->{'annot'}->{$append_file} ) ) { 
                    $tab_count = $index_lines_ref->{'annot_file'}->{$append_file}->{'tab_count'};
                    my $added_tabs = ("\t" x $tab_count);
                    print "$added_tabs";
                }
            }
            print "\n";
        }
        if ($semicolon) { 
            if ($full_annot) { 
                die "Not yet able to handle both semicolons and complex annotations!\n";
            }
            my $text = $index_lines_ref->{'key_id'}->{$annot_key_id}->{'text'};
            my $endtext = q{};
            if ( $text =~ / \A (.*?) ([.]{0,1} ["]{0,1}) \s* \z /xms ) { 
                $text    = $1;
                $endtext = $2;
                print "$text";
                foreach my $append_file (@annot_files) { 
                    if (exists $index_lines_ref->{'key_id'}->{$annot_key_id}->{'annot'}->{$append_file} ) { 
                        print "; $index_lines_ref->{'key_id'}->{$annot_key_id}->{'annot'}->{$append_file}";
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

