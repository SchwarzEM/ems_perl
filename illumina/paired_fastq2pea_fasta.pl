#!/usr/bin/env perl

# paired_fastq2pea_fasta.pl -- Erich Schwarz <emsch@its.caltech.edu>, 1/13/2011.
# Purpose:  From 1+ files/input stream, take eight lines at a time.  Check for superficially correct headers and enforce paired-end-ness.  Convert from FASTQ to PE-Assembler-style "fasta" file, in which second read is reverse-complemented and listed after the first read and a "\t" spacer.
# Documentation for correct output format:  http://www.comp.nus.edu.sg/~bioinfo/peasm/PE_manual.htm; also, undocumented requirement for all lowercase.
# General documentation: http://bioinformatics.oxfordjournals.org/content/early/2010/12/12/bioinformatics.btq626.abstract

use strict;
use warnings;
use Getopt::Long;

my @input_files = ();
my $input       = q{};
my %opts        = ();
my $read_name   = q{};
my $read_stem   = q{};

my $data_ref;

# Default suffixes for paired-end read 1 and paired-end read 2.
# Note that these are written '\/' because, 
#      although they are literal strings, they'll gointo a 
#      pattern-match which requires that '/' be an escape character.

$opts{'r1'} = '#0\/1';
$opts{'r2'} = '#0\/2';

GetOptions(
    'input=s{,}'        => \@input_files,
    'r1=s'              => \$opts{'r1'},
    'r2=s'              => \$opts{'r2'}, 
    'help'              => \$opts{'help'},
);

if ( $opts{'help'} or (! @input_files) or (! $opts{'r1'} ) or (! $opts{'r2'} ) ) {
    print "\n";
    print "usage: paired_fastq2pea_fasta.pl\n";
    print "       --input|-i        [input file(s); accepts '-' for input stream]\n";
    print "       --r1 <suffix1>    [suffix marking paired-end read 1, default \"#0/1\"]\n";
    print "       --r2 <suffix2>    [suffix marking paired-end read 2, default \"#0/1\"]\n";
    print "       --help|-h         [help - print this message]\n";
    print "\n";
    exit;
}

# Define reusable filehandle, and initialize line-counting; 
#     line-counting to be continued through all input files.

my $INPUT_FILE;
my $i = 0;

foreach my $input_file (@input_files) { 
    # Accept input data either from a stream (via '-') or from a standard file.
    if ($input_file eq '-') {
        # Special case: get the stdin handle
        $INPUT_FILE = *STDIN{IO};
    }
    else {
        # Standard case: open the file
        open $INPUT_FILE, '<', $input_file or die "Can't open $input_file. $!\n";
    }

    # Process data, keeping steady track of lines in all successive input files.
    while (my $input = <$INPUT_FILE>) {
        chomp $input;
        $i++;
        my $j = ( $i % 8);
        if ( $j == 1 ) { 

            # Avoid normal 'xms' parsing to prevent ambiguous regexes.
            if ( $input !~ /\A[@]\S+/ ) { 
                die "Can't parse putative header: $input\n";
            }

            $read_name = q{};
            if ( $input =~ /\A[@](\S+)/ ) {
                $read_name = $1;
            }

            if ( $read_name !~ /\A\S+$opts{'r1'}\z/ ) { 
                die "Can't parse read name $read_name in header: $input\n";
            }
            elsif ( ( $read_name =~ /\A\S+$opts{'r1'}\z/ ) and ( $read_name =~ /\A\S+$opts{'r2'}\z/ ) ) {
                die "Header $input somehow has both suffixes $opts{'r1'} and $opts{'r2'}!\n";
            }
            elsif ( $read_name =~ /\A(\S+)$opts{'r1'}\z/ ) {
                # Get the stem name for an r1 read.
                $read_stem = $1;

                # Die loudly if an identically named r1 has already been read.
                if ( exists $data_ref->{'observed_name'}->{$read_stem}->{'r1'} ) {
                    die "Multiple sequence entries with sequence name $read_stem$opts{'r1'}\n";
                }

                # Die loudly if its r2 already has been read.
                if ( exists $data_ref->{'observed_name'}->{$read_stem}->{'r2'} ) {
                    die "Reads for $read_stem are in wrong order!\n";
                }

                # Record the reading of first (and only allowed) r1 entry.
                $data_ref->{'observed_name'}->{$read_stem}->{'r1'} = 1;
            }
            else { 
                die "Can't parse input: $input\n";
            }
        }

        # Get the FASTA sequence associated with namestem/r1 read.
        elsif ( $j == 2 ) { 
            $input =~ s/\s+\z//;
            if ( $input !~ /\A[acgtnACGTN]+\z/ ) {
                die "Can't parse sequence of $read_name: $input\n";
            }
            $data_ref->{'observed_name'}->{$read_stem}->{'r1'} = $input;
        }

        # On line 5 of 8, require an r2 match for the namestem.
        elsif ( $j == 5 ) { 

            my $curr_read_stem = q{};
            if ( $input !~ /\A[@]\S+/ ) {
                die "Can't parse putative header: $input\n";
            }
    
            $read_name = q{};
            if ( $input =~ /\A[@](\S+)/ ) {
                $read_name = $1;
            }

            if ( $read_name !~ /\A\S+$opts{'r2'}\z/ ) {
                die "Can't parse read name $read_name in header: $input\n";
            } 
            elsif ( ( $read_name =~ /\A\S+$opts{'r1'}\z/ ) and ( $read_name =~ /\A\S+$opts{'r2'}\z/ ) ) {
                die "Header $input somehow has both suffixes $opts{'r1'} and $opts{'r2'}!\n";
            }
            elsif ( $read_name =~ /\A(\S+)$opts{'r2'}\z/ ) {
                # Get the stem name for an r1 read.
                $curr_read_stem = $1;

                # Enforce identical read stems between pairs:
                if ( $curr_read_stem ne $read_stem ) { 
                    die "Read stems $read_stem and $curr_read_stem do not match!\n";
                }

                # Die loudly if an identically named r2 has already been read.
                if ( exists $data_ref->{'observed_name'}->{$read_stem}->{'r2'} ) {
                    die "Multiple sequence entries with sequence name $read_stem$opts{'r1'}\n";
                }

                # Die loudly if its r1 has not been read.
                if (! exists $data_ref->{'observed_name'}->{$read_stem}->{'r1'} ) {
                    die "Missing first read for $read_stem!\n";
                }

                # Record the reading of first (and only allowed) r2 entry.
                $data_ref->{'observed_name'}->{$read_stem}->{'r2'} = 1;
            }
            else {
                die "Can't parse input: $input\n";
            }
        }

        # Get the FASTA sequence associated with namestem/r12 read, and print out the appropriate pseudo-FASTA.
        elsif ( $j == 6 ) { 
            $input =~ s/\s+\z//;
            if ( $input !~ /\A[acgtnACGTN]+\z/ ) {
                die "Can't parse sequence of $read_name: $input\n";
            }
            if ( $input =~ /\A[acgtnACGTN]+\z/ ) {
                my $second_read   = revcomp($input);
                my $header_line   = '>' . "$read_stem";
                my $sequence_line = "$data_ref->{'observed_name'}->{$read_stem}->{'r1'}" 
                                    . "\t"
                                    . "$second_read"
                                    ;
                $sequence_line =~ tr/[A-Z]/[a-z]/;
                if ( $sequence_line =~ / (n | [.]) /xms ) { 
                    die "PE-Assembler cannot accept 'n' or '.' character in sequence line:\n $sequence_line\n";
                }
                print "$header_line\n";
                print "$sequence_line\n";
            }
            $read_name = q{};  
            $read_stem = q{};
        }
    }
    close $INPUT_FILE or die "Can't close filehandle to $input_file: $!\n";
}

# Note that -- appropriately -- this does nothing to 'n' or 'N' residues.
sub revcomp { 
    my $in_string = $_[0];
    $in_string =~ tr/[acgtACGT]/[tgcaTGCA]/;
    my @in_residues = split //, $in_string;
    @in_residues = reverse @in_residues;
    my $out_string = join q{}, @in_residues;
    return $out_string;
}

