#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

my $data_ref;
my $header = q{};

my @infiles = ();
my $qthresh = 1;
my $help;

GetOptions ( 'infiles=s{,}' => \@infiles,
             'qthreshold=f' => \$qthresh,
             'help'         => \$help,   );

if ( $help or (! @infiles) or (! looks_like_number($qthresh) ) or ( $qthresh < 0 ) or ( $qthresh > 1 ) ) { 
    die "Format: make_edgeR_result_summary_28apr2013.pl\n",
        "    --infile|-i      <input stream/files>\n",
        "    --qthreshold|-q  [q-value threshold between 0 and 1; default is 1]\n",
        "    --help|-h        [print this message]\n",
        ;
}

foreach my $infile (@infiles) {
    my $INPUT_FILE;
    if ($infile eq '-') {
        # Special case: get the stdin handle
        $INPUT_FILE = *STDIN{IO};
    }
    else {
        # Standard case: open the file
        open $INPUT_FILE, '<', $infile or die "Can't open input file $infile. $!\n";
    }
    while (my $input = <$INPUT_FILE>) { 

# 4sp.

        chomp $input;

        # This conditional will only be invoked once, when reading the first line of the input table.
        # Either that first line gets parsed into satisfactory fields, or the script dies.
        if (! exists $data_ref->{'column'} ) { 
            my @input_fields = split /\t/, $input;

            # Reality-check the header's fields.

            # Capture original column number:
            my $input_field_count = @input_fields;

            # Test for Gene in first column, and, simultaneously, remove it from @input_fields:
            my $first_header_value = shift @input_fields;
            if ( $first_header_value ne 'Gene') {
                die "Can't find \"Gene\" in first column of header: $input\n";
            }

            # Test for Gene + even no. columns.
            $input_field_count--;
            if ( ( $input_field_count % 2 ) != 0 ) {
                die "Odd number of data columns after \"Gene\" in header: $input\n";
            }

            # Step through the columns one by one.  Enforce pairs like this:
            # ACEY.L3i.vs.ACEY.24.PI.edgeR.logFC      ACEY.L3i.vs.ACEY.24.PI.edgeR.q-value
            my $i          = 0;
            my $comparison = q{};
            foreach my $input_field (@input_fields) {
                $i++;
                if (! $comparison ) { 
                    if ( $input_field =~ /\A (\S+ \.vs\. \S+) \.edgeR\.logFC \z/xms ) { 
                        $comparison = $1;
                        $data_ref->{'column'}->{$i}->{'comparison'} = $comparison;
                        $data_ref->{'column'}->{$i}->{'content'}    = 'fold_change';
                    }
                    else { 
                        die "Can't parse input field number $i from header: $input\n";
                    }
                }
                else { 
                    if ( $input_field =~ /\A (\S+ \.vs\. \S+) \.edgeR\.q-value \z/xms ) {
                        my $comparison_new = $1;
                        if ( $comparison ne $comparison_new ) { 
                            die "Inconsistent comparisons ($comparison vs. $comparison_new) in input: $input\n";
                        }
                        $data_ref->{'column'}->{$i}->{'comparison'} = $comparison;
                        $data_ref->{'column'}->{$i}->{'content'}    = 'q-value';
                        $comparison = q{};
                    }
                    else {
                        die "Can't parse input field number $i from header: $input\n";
                    }
                }
            }
        }

        # This 'else' can only occur if (exists $data_ref->{'column'}).
        # After reading the first line of the file, the script should always default to here.

        else {
            my @input_fields = split /\t/, $input;

            my $gene = shift @input_fields;
            if ( $gene !~ /\S/xms ) { 
                die "Need at least one non-space character in gene \"$gene\", in input: $input\n";
            }
            if ( exists $data_ref->{'gene'}->{$gene} ) { 
                die "Redundant gene \"$gene\", in input: $input\n";
            }

            my $i = 0;
            foreach my $input_field (@input_fields) { 
                $i++;
                if (! looks_like_number($input_field) ) { 
                    die "Non-numerical input field \"$input_field\" in input: $input\n";
                }
                # E.g., $comparison = 'ACEY.L3i.vs.ACEY.24.PI':
                my $comparison = $data_ref->{'column'}->{$i}->{'comparison'};

                # Will be either 'fold_change' or 'q-value':
                my $content = $data_ref->{'column'}->{$i}->{'content'};

                # Note: we *will* get some cases of 0 fold change.  
                # These will be q_value of 1, of course, but if we allow q = 1, the edge case must be dealt with.
                # In practice we'll give them a class of their own -- '0' rather than '+' or '-'!

                # Record the data:
                $data_ref->{'gene'}->{$gene}->{$comparison}->{$content} = $input_field;
            }
        }
    }
    close $INPUT_FILE or die "Can't close filehandle to input file $infile. $!\n";
}


my @genes = sort keys %{ $data_ref->{'gene'} };
if (@genes) { 
    $header = "Gene\tedgeR";
}

foreach my $gene (@genes) {
    my @instances = ();

    # Sort comparisons from q = 0 to q > 0. 
    my @comparisons = grep { $data_ref->{'gene'}->{$gene}->{$_}->{'q-value'} <= $qthresh }
                      sort {     $data_ref->{'gene'}->{$gene}->{$a}->{'q-value'} 
                             <=> $data_ref->{'gene'}->{$gene}->{$b}->{'q-value'} } 
                      keys %{ $data_ref->{'gene'}->{$gene} };
    foreach my $comparison (@comparisons) {
        my $fold_change = $data_ref->{'gene'}->{$gene}->{$comparison}->{'fold_change'};
        my $fold_sign = q{};

        if ( $fold_change > 0 ) {
            $fold_sign = '+';
        }
        if ( $fold_change < 0 ) { 
            $fold_sign = '-';
        }
        # Only relevant if we are really not trying to filter out q = 1 genes...
        if ( $fold_change == 0 ) {
            $fold_sign = '0';
        }

        if ( $comparison =~ /\A (\S+) \.vs\. (\S+) \z/xms ) {
            my $condition_a = $1;
            my $condition_b = $2;
            my $q_value     = $data_ref->{'gene'}->{$gene}->{$comparison}->{'q-value'};
            my $instance    = "$condition_a to $condition_b [$fold_sign: $q_value]";
            push @instances, $instance;
        }
    }

    my $instance_text = q{};
    if (@instances) {
        $instance_text = join '; ', @instances;
    }

    # Print header once:
    print "$header\n" if $header;
    $header = q{};

    print "$gene\t$instance_text\n";
}

