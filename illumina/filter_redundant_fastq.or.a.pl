#!/usr/bin/env perl

# filter_redundant_fastq.or.a.pl -- Erich Schwarz <emsch@caltech.edu>, 9/2/2011.
# Purpose: given a stream of Illumina data in either FASTQ or FASTA format, pass only reads with non-redundant names (with no effort to check contents).

my $i     = 0;
my @lines = ();
my $head1 = q{};
my $text1 = q{};
my $head2 = q{};
my $text2 = q{};
my $read  = q{};
my %seen  = ();

while (my $input = <>) { 
    # Count lines by human 1-based numbering, not Perl 0-based:
    $i++;
    
    # Each stanza of four lines resets the line storage array:
    if ( ( $i % 4 ) == 1 ) { 
        @lines = ();
        # No carry-over of other variables:
        $head1 = q{};
        $text1 = q{};
        $head2 = q{};
        $text2 = q{};
        $read  = q{};
    }

    # Regardless of the state of the line-storage array, do this for every single input line:
    chomp $input;
    push @lines, $input;

    # When modulus reaches 0 (i.e., a multiple of 4), it is time to process text:
    if ( ( $i % 4 ) == 0 ) {
        # Less efficient than just using $lines[0] etc., but more readable:
        $head1 = $lines[0];
        $text1 = $lines[1];
        $head2 = $lines[2];
        $text2 = $lines[3];

        # Default to assuming FASTQ input if two very basic requirements are met.
        if ( ( $text1 =~ /\A [ACGTNacgtn]+ \z/xms ) and ( $head1 =~ /\A [@] (\S+) /xms ) ) { 
            $read = $1;
            if (! exists $seen{$read} ) { 
                print "$head1\n";
                print "$text1\n";
                print "$head2\n";
                print "$text2\n";
                $seen{$read} = 1;
            }
        }

        # Opt instead for FASTA input if lines clearly fit FASTA instead.
        elsif (     ( $head1 =~ /\A [>] \S+ /xms )
                and ( $text1 =~ /\A [ACGTNacgtn]+ \z/xms ) 
                and ( $head2 =~ /\A [>] \S+ /xms )
                and ( $text2 =~ /\A [ACGTNacgtn]+ \z/xms ) ) { 
            foreach my $data_ref ( [$head1, $text1], [$head2, $text2] ) { 
                if ( ( $data_ref->[1] =~ /\A [ACGTNacgtn]+ \z/xms ) and ( $data_ref->[0] =~ /\A [>] (\S+) /xms ) ) { 
                    $read = $1;
                    if (! exists $seen{$read} ) {
                        print "$data_ref->[0]\n";
                        print "$data_ref->[1]\n";
                        $seen{$read} = 1;
                    }
                }
            }
        }

        # If parsing fails, the program dies loudly:
        else { 
            die "Can't parse lines:\n",
                "\n",
                "$head1\n",
                "$text1\n",
                "$head2\n",
                "$text2\n",
                ;
        }
    }
}

