#!/usr/bin/env perl

# make_max_obs_data_table.pl -- Erich Schwarz <ems394@cornell.edu>, 2/14/2013.
# Purpose: given numerical data tables (e.g., 2+ ".results" tables from RSEM 1.2.0 for genes or isoforms), compile their maximum observed values for each gene.

use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my @infiles        = ();
my $general_header = q{};
my $id_type        = ();
my @data_fields    = ();
my $data_tabcount  = 0;
my $data_ref;
my $rsem;
my $help;

GetOptions ( 'infiles=s{,}' => \@infiles,
             'rsem'         => \$rsem,
             'help'         => \$help,   );

if ( $help or (! @infiles) ) { 
    die "Format: make_max_obs_data_table.pl\n",
        "    --infile|-i   [input files]\n",
        "    --rsem|-r     [deal with RSEM input, which has two left-hand id columns (main and subsidiary)]\n",
        "    --help|-h     [print this message]\n",
        ;
}

foreach my $infile (@infiles) {
    open my $INPUT_FILE, '<', $infile or die "Can't open input file $infile. $!\n";
    my $file_header = q{};
    while (my $input = <$INPUT_FILE>) { 
        chomp $input;

        # This will only happen on the first line of a new file, i.e., before we've recorded its header.
        if (! $file_header ) {
            # Define that file's header.
            $file_header = $input;

            # This will only happen on the first line of the first file, i.e., before we've defined its header as the required standard.
            if (! $general_header) { 
                $general_header = $file_header;

                # When we define the general header text, enforce its format, name its data fields, and count its tabs.
                if ( $general_header =~ /\A (\S+) \t ( (?: \S+ \t )* \S+ ) \z /xms ) {
                    $id_type       = $1;
                    my $data_text  = $2;
                    @data_fields   = split "\t", $data_text;
                    $data_tabcount = ( $data_text =~ tr/\t/\t/ );
                }
                else { 
                    die "From infile $infile, can't parse header format: $general_header\n";
                }
            }
            # If we had previously recorded a header, require identity to each new file's header.
            elsif ( $file_header ne $general_header ) {
                die "From infile $infile, inconsistent header format:\n", "$general_header\n", "vs.\n", "$file_header\n", ;
            }
        }
        else { 
            if ( $input =~ /\A (\S+) \t ( (?: \S+ \t ){$data_tabcount} \S+ ) \z /xms ) {
                my $id         = $1;
                my $data_text  = $2;
                my @data       = split "\t", $data_text;
                foreach my $i (0..$data_tabcount) { 
                    my $datum = $data[$i];
                    my $type  = $data_fields[$i];

                    # For data that we know are supposed to be numerical, enforce number and pick largest value seen so far.
                    if ( ( $i >= 1 ) or (! $rsem ) ) {
                        if (! looks_like_number($datum) ) { 
                            die "From infile $infile, data point $datum not numeric in data line: $input\n",
                                "For RSEM tables, try --rsem|-r\n",
                                ;
                        } 
                        if ( (! exists $data_ref->{'id'}->{$id}->{'data_type'}->{$type} ) or ( $datum > $data_ref->{'id'}->{$id}->{'data_type'}->{$type} ) ) { 
                            $data_ref->{'id'}->{$id}->{'data_type'}->{$type} = $datum;
                        }
                    }
                    # For non-numerical data of a gene, we assume that all tables have the same value, allowing us to update that value for every data file.
                    else { 
                        $data_ref->{'id'}->{$id}->{'data_type'}->{$type} = $datum;
                    }
                }
            }
            else { 
                die "From infile $infile, can't parse data line: $input\n";
            }
        }
    }
    close $INPUT_FILE or die "Can't close filehandle to input file $infile. $!\n";
}

print "$general_header\n";

my @id_set = sort keys %{ $data_ref->{'id'} };
foreach my $id (@id_set) { 
    my @output_fields = ();
    push @output_fields, $id;
    foreach my $type (@data_fields) { 
        my $datum = $data_ref->{'id'}->{$id}->{'data_type'}->{$type};
        push @output_fields, $datum;
    }
    my $output_text = join "\t", @output_fields;
    print "$output_text\n";
}

