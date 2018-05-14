#!/usr/bin/env perl

# paired_vs_unp_fastq.or.a_legacy.pl -- Erich Schwarz <emsch@its.caltech.edu>, 2/28/2011, ARCHIVAL LEGACY version -- kept so that I can (more or less) reproduce past work; note that I fixed one visible bug in this version that had been present in past work.
# Purpose:  Take two or four lines at a time, depending on FASTA or (default) FASTQ input.  Check for superficially correct headers.  Namesort (avoiding RAM and brute force) into paired and unpaired output files.

use strict;
use warnings;
use Getopt::Long;

my $input = q{};

my $head1  = q{};
my $nt_seq = q{};
my $head2  = q{};
my $quals  = q{};

my %opts  = ();

my $data_ref;
my $read_stem = q{};

my $fasta_input = 0;
my $fastq_input = 0;

# Default suffixes for paired-end read 1 and paired-end read 2:
$opts{'r1'} = '#0\/1';
$opts{'r2'} = '#0\/2';

GetOptions(
    'input=s'           => \$opts{'input'},
    'paired_output:s'   => \$opts{'paired'},  
    'unpaired_output:s' => \$opts{'unpaired'},
    'fasta'             => \$fasta_input,
    'fastq'             => \$fastq_input,
    'r1=s'              => \$opts{'r1'},
    'r2=s'              => \$opts{'r2'}, 
    'help'              => \$opts{'help'},
);

if (    (! $opts{'input'}) 
     or (! $opts{'r1'}   ) 
     or (! $opts{'r2'}   ) 
     or ( $fasta_input and $fastq_input ) ) {
    $opts{'help'} = 1;
}

if ( $opts{'help'} ) {
    print "\n";
    print "usage: paired_vs_unp_fastq.or.a.pl\n";
    print "       -i|--input            <name1>    input file name (fasta or fastq), mandatory\n";
    print "       -p|--paired_output    <name2>    output paired-end file name (fastq); default \"name1.paired\"\n";
    print "       -u|--unpaired_output  <name3>    output unpaired-end file name (fastq); default \"name1.unpaired\"\n";
    print "       --fastq                          fastQ format for input and outputs (default)\n";
    print "       --fasta                          fastA format for input and outputs (mutually exclusive with fastQ)\n";
    print "       --r1                  <suffix1>  suffix marking paired-end read 1, default \"#0/1\".\n";
    print "       --r2                  <suffix2>  suffix marking paired-end read 2, default \"#0/2\".\n";
    print "       -h|--help                        help - print this message\n";
    print "\n";
    exit;
}

if ( (! $fasta_input ) and (! $fastq_input ) ) { 
    $fastq_input = 1;
}

if (! $opts{'paired'} ) { 
    $opts{'paired'} = $opts{'input'} . ".paired";
    $opts{'paired'} = safename( $opts{'paired'} );
}

if (! $opts{'unpaired'} ) {
    $opts{'unpaired'} = $opts{'input'} . ".unpaired";
    $opts{'unpaired'} = safename( $opts{'unpaired'} );
}

# Do this first, to ensure that the output files *can* be opened (otherwise, die preemptively):
open my $PAIRED,   '>', $opts{'paired'}   or die "Cannot open output paired-reads file $opts{'paired'}: $!";  
open my $UNPAIRED, '>', $opts{'unpaired'} or die "Cannot open output unpaired-red file $opts{'unpaired'}: $!";

# Then, go through input file completely to get accurate classification of reads (before even trying to deal with contents).
open my $INFILE,   '<', $opts{'input'}    or die "Cannot open input jumbled-reads file $opts{'input'}: $!";

my $i = 0;
my $j = 0;
my $lead_pattern = q{};
if ( $fastq_input ) {
    $lead_pattern = '[@]';
}
if ( $fasta_input ) {
    $lead_pattern = '>';
}

while (my $input = <$INFILE>) {
    chomp $input;
    $i++;
    if ( $fastq_input ) {
        $j = ( $i % 4);
    }    
    if ( $fasta_input ) {
        $j = ( $i % 2);
    }    
    if ($j == 1) { 

        # Avoid normal 'xms' parsing to prevent ambiguous regexes.

        if ( $input !~ /\A$lead_pattern\S+/ ) { 
            die "Can't parse putative header: $input\n";
        }

        my $read_name = q{};
        if ( $input =~ /\A$lead_pattern(\S+)/ ) {
            $read_name = $1;
        }

        if ( ( $read_name !~ /\A\S+$opts{'r1'}\z/ ) and ( $read_name !~ /\A\S+$opts{'r2'}\z/ ) ) { 
            die "Can't parse read name $read_name in header: $input\n";
        }
        elsif ( ( $read_name =~ /\A\S+$opts{'r1'}\z/ ) and ( $read_name =~ /\A\S+$opts{'r2'}\z/ ) ) {
            die "Header $input somehow has both suffixes $opts{'r1'} and $opts{'r2'}!\n";
        }
        else { 
            if ( $read_name =~ /\A(\S+)$opts{'r1'}\z/ ) {
                # Get the stem name for an r1 read.
                $read_stem = $1;

                # Die loudly if an identically named r1 has already been read.
                if ( exists $data_ref->{'observed_name'}->{$read_stem}->{'r1'} ) {
                    die "Multiple sequence entries with sequence name $read_stem$opts{'r1'}\n";
                }

                # Record the reading of first (and only allowed) r1 entry.
                $data_ref->{'observed_name'}->{$read_stem}->{'r1'} = 1;

                # Classify as jumbled, not singleton, if its r2 already has been read.
                if ( exists $data_ref->{'observed_name'}->{$read_stem}->{'r2'} ) {
                    $data_ref->{'jumbled'}->{$read_stem} = 1;

                    # If it's jumbled, it's no longer a singleton, and it can't be paired with the preceding read:
                    delete $data_ref->{'singleton'}->{$read_stem} ;
                    delete $data_ref->{'maybe_paired_properly'};
                }

                # Classify as singleton, and as *possibly* correctly paired, if its r2 has not yet been read.
                if (! exists $data_ref->{'observed_name'}->{$read_stem}->{'r2'} ) { 
                    $data_ref->{'singleton'}->{$read_stem} = 1;

                    # The following step automatically wipes out the *last* read stem listed as 'maybe_paired_properly'.
                    # Note that keeping this straight demands remembering the *immediate* predecessor, not anything earlier.
                    $data_ref->{'maybe_paired_properly'} = $read_stem;
                }
            }
            elsif ( $read_name =~ /\A(\S+)$opts{'r2'}\z/ ) {
                # Get the stem name for an r2 read.
                $read_stem = $1;

                # Die loudly if an identically named r2 has already been read.
                if ( exists $data_ref->{'observed_name'}->{$read_stem}->{'r2'} ) {
                    die "Multiple sequence entries with sequence name $read_stem$opts{'r2'}\n";
                }

                # Note: while editing the above in a successor script, I found that this conditional was actually:
                # "if ( exists $data_ref->{$read_stem}->{'r2'} ) {"
                # which is just wrong!  Hopefully this error hasn't led to too many actual foul-ups of filtered FASTQ files...
                # At any rate, I've corrected the error in this archival legacy version.

                # Record the reading of first (and only allowed) r2 entry.
                $data_ref->{'observed_name'}->{$read_stem}->{'r2'} = 1;

                # Classify as singleton, but as unpairable with previous read, if its r1 has not yet been read.
                if (! exists $data_ref->{'observed_name'}->{$read_stem}->{'r1'} ) {
                    $data_ref->{'singleton'}->{$read_stem} = 1;
                    delete $data_ref->{'maybe_paired_properly'};
                }

                # Classify the r2 as coming after its r1.
                if ( exists $data_ref->{'observed_name'}->{$read_stem}->{'r1'} ) { 
                    # *If* the read just before this one was correct ...
                    if ( $data_ref->{'maybe_paired_properly'} and ( $read_stem eq $data_ref->{'maybe_paired_properly'} ) ) { 
                        # ... classify as correctly paired!
                        $data_ref->{'paired'}->{$read_stem} = 1;

                        # And immediately wipe out these annotations:
                        delete $data_ref->{'singleton'}->{$read_stem} ;
                        delete $data_ref->{'maybe_paired_properly'};
                    }
                    # Otherwise, it is jumbled, not singletone or pairable.
                    else { 
                        $data_ref->{'jumbled'}->{$read_stem} = 1;
                        delete $data_ref->{'singleton'}->{$read_stem} ;
                        delete $data_ref->{'maybe_paired_properly'};
                    }
                }
            }
            # Compel deterministic parsing.
            else { 
                die "Cannot decipher header: $input\n";
            }
        }
    }
}
close $INFILE, or die "Cannot open close filehandle to jumbled-reads file $opts{'input'}: $!";

# Finally, go through input file and sort out data.
open $INFILE, '<', $opts{'input'}    or die "Cannot open input jumbled-reads file $opts{'input'}: $!";
while (<$INFILE>) { 
    # Tricky Perl syntax: 
    #     it forces me to start with $_ already loaded from <$INFILE>,
    #     but afterwards needs prompting for linefeeds.

    $head1 = $_;
    chomp $head1;
    $nt_seq = <$INFILE>;
    chomp $nt_seq;
    if ( $fastq_input ) { 
        $head2  = <$INFILE>;
        chomp $head2;
        $quals  = <$INFILE>;
        chomp $quals;
    }

    # Overspecify this to avoid confusing the Perl interpreter.
    my @read_data     = ( $head1, $nt_seq, );
    if ( $fastq_input ) {
        push @read_data, ( $head2, $quals, );
    }
    my $read_data_ref = \@read_data;

    my $read_name = q{};
    if ( $head1 =~ /\A$lead_pattern(\S+)/ ) { 
        $read_name = $1;
    }
    if ( $read_name =~ /\A(\S+)(?:$opts{'r1'}|$opts{'r2'})\z/ ) {
        $read_stem = $1;
        if ( exists $data_ref->{'singleton'}->{$read_stem} ) { 
            my $print_text = join "\n", @read_data;
            print $UNPAIRED "$print_text\n";
        }
        elsif ( exists $data_ref->{'paired'}->{$read_stem} ) { 
            my $print_text = join "\n", @read_data;
            print $PAIRED "$print_text\n";
        }
        elsif ( exists $data_ref->{'jumbled'}->{$read_stem} ) { 
            if ( $head1 =~ /\A$lead_pattern(\S+)$opts{'r1'}\z/ ) {
                $data_ref->{'data'}->{$read_stem}->{'r1'} = $read_data_ref;
            }
            elsif ( $head1 =~ /\A$lead_pattern(\S+)$opts{'r2'}\z/ ) {
                $data_ref->{'data'}->{$read_stem}->{'r2'} = $read_data_ref;
            }
            else {
                die "Unable to resolve which paired-end the jumbled read $head1 was.\n";
            }
        }
        else { 
            die "Unable to resolve status of read $head1.\n";
        }
    }
    else { 
        die "From input file $opts{'input'}, unable to parse: $head1\n";
    }
}
close $INFILE,  or die "Cannot open close filehandle to jumbled-reads file $opts{'input'}: $!";
close $UNPAIRED or die "Cannot close filehandle to unpaired-reads file $opts{'unpaired'}: $!";

my @read_stems = sort keys %{ $data_ref->{'data'} };
foreach my $read_stem2 (@read_stems) { 
    # Set to print to the paired file, if both ends of a read are available:
    if ( ( exists $data_ref->{'data'}->{$read_stem2}->{'r1'} ) and ( exists $data_ref->{'data'}->{$read_stem2}->{'r2'} ) ) {
        my $print_text = join "\n", @{ $data_ref->{'data'}->{$read_stem2}->{'r1'} };
        print $PAIRED "$print_text\n";
        $print_text = join "\n", @{ $data_ref->{'data'}->{$read_stem2}->{'r2'} };
        print $PAIRED "$print_text\n";
    }
    else { 
        die "Unable to properly parse jumbled reads stored in memory for sorted output to $opts{'paired'}\n";
    }
}
close $PAIRED or die "Cannot close filehandle to paired-reads file $opts{'paired'}: $!";


# Purpose: always print/export data to a filename that's new, avoiding overwriting existing files.
sub safename {
    my $filename = $_[0];
    my $orig_filename = $filename;
    if (-e $orig_filename) {
        my $suffix1 = 1;
        $filename = $filename . ".$suffix1";
        while (-e $filename) {
            $suffix1++;
            $filename =~ s/\.\d+\z//xms;
            $filename = $filename . ".$suffix1";
        }
    }
    return $filename;
}

