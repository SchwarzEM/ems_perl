#!/usr/bin/env perl

# split_paired_fastq.or.a.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/9/2011.
# Purpose:  Depending on input (FASTA or FASTQ), take two or four lines at a time.  Check for superficially correct headers.  Split into two end-1 and end-2 files.

use strict;
use warnings;
use Getopt::Long;

my $input = q{};

my $head1  = q{};
my $nt_seq = q{};
my $head2  = q{};
my $quals  = q{};

my %opts  = ();

my $fasta_input = 0;
my $fastq_input = 0;

my $read_stem = q{};

# Default suffixes for paired-end read 1 and paired-end read 2:
$opts{'r1'} = '#0\/1';
$opts{'r2'} = '#0\/2';

GetOptions(
    "input=s" => \$opts{'input'},
    "o1:s"    => \$opts{'output1'},  
    "o2:s"    => \$opts{'output2'},
    "r1=s"    => \$opts{'r1'},
    "r2=s"    => \$opts{'r2'}, 
    "q|fastq" => \$fastq_input,
    "a|fasta" => \$fasta_input,
    "help"    => \$opts{'help'},
);

if (    (! $opts{'input'}               ) 
     or (! $opts{'r1'}                  ) 
     or (! $opts{'r2'}                  ) 
     or ( $fastq_input and $fasta_input ) ) {
    $opts{'help'} = 1;
}

if ( $opts{'help'} ) {
    print "usage: split_paired_fastq.or.a.pl\n";
    print "       -i|--input  <name1>    input file name (fastq/a), mandatory (though '-' will work)\n";
    print "       -q|--fastq             fastQ format for input and outputs (default)\n";
    print "       -a|--fasta             fastA format for input and outputs (mutually exclusive with fastQ)\n";
    print "       --o1        <name2>    output paired-end 1 file name (fastq/a); default \"name1.output1\"\n";
    print "       --o2        <name3>    output paired-end 2 file name (fastq/a); default \"name1.output2\"\n";
    print "       --r1        <suffix1>  suffix marking paired-end read 1, default \"#0\\/1\".\n";
    print "       --r2        <suffix2>  suffix marking paired-end read 2, default \"#0\\/2\".\n";
    print "       -h|--help              help - print this message\n";
    exit;
}

if ( (! $fasta_input ) and (! $fastq_input ) ) { 
    $fastq_input = 1;
}

if (! $opts{'output1'} ) { 
    $opts{'output1'} = $opts{'input'} . '.pe1';
    $opts{'output1'} = safename( $opts{'output1'} );
}

if (! $opts{'output2'} ) {
    $opts{'output2'} = $opts{'input'} . '.pe2';
    $opts{'output2'} = safename( $opts{'output2'} );
}

# Accept either a stream from '-' or a standard input file.
my $INFILE;
if ($opts{'input'} eq '-') {
    # Special case: get the stdin handle
    $INFILE = *STDIN{IO};
}
else {
    # Standard case: open the file
    open $INFILE, '<', $opts{'input'}   or die "Cannot open input jumbled-reads file $opts{'input'}: $!"; 
}

# Do this first, to ensure that outputs can happen at all.
open my $OUTFILE_1, '>', $opts{'output1'} or die "Cannot open output paired-end-1 file $opts{'output1'}: $!";
open my $OUTFILE_2, '>', $opts{'output2'} or die "Cannot open output paired-end-2 file $opts{'output2'}: $!";

my $lead_pattern = q{};
if ( $fastq_input ) {
    $lead_pattern = '[@]';
}
if ( $fasta_input ) {
    $lead_pattern = '>';
}

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

        # Note: for reasons that are highly unclear to me as of 10/18/2010, the normally-good
        #    'xms' modifier on regular expressions causes parsing of header names from their
        #    paired-end suffixes to totally fail!
        #    So, to get the script to simply work, primitive non-legible non-xms regexes are
        #        used whenever FASTQ headers and suffixes need to be sorted out.
    
        if ( $head2 !~ /\A \+ /xms ) { 
            die "Can't parse header 2: $head2\n";
        }
    }

    if (     ( $head1 !~ /\A$lead_pattern\S+$opts{'r1'}\z/ )
         and ( $head1 !~ /\A$lead_pattern\S+$opts{'r2'}\z/ ) ) { 
        die "Can't parse header 1: $head1\n";
    }
    elsif (     ( $head1 !~ /\A$lead_pattern\S+$opts{'r1'}\z/ ) 
            and ( $head1 !~ /\A$lead_pattern\S+$opts{'r2'}\z/ ) ) { 
        die "Header $head1 somehow has both suffixes $opts{'r1'} and $opts{'r2'}!\n";
    }
    else { 
        my @output_lines = ( $head1, $nt_seq, );
        if ( $fastq_input ) { 
            push @output_lines, ( $head2, $quals, );
        }        
        my $print_text = join "\n", @output_lines;

        if ( $head1 =~ /\A$lead_pattern(\S+)$opts{'r1'}\z/ ) {
            print $OUTFILE_1 "$print_text\n";
        }
        elsif ( $head1 =~ /\A$lead_pattern(\S+)$opts{'r2'}\z/ ) {
            print $OUTFILE_2 "$print_text\n";
        } 
    }
}
close $INFILE, or die "Cannot open close filehandle to jumbled-reads file $opts{'input'}: $!";
close $OUTFILE_1 or die "Cannot close filehandle to paired-end-1 file $opts{'output1'}: $!";
close $OUTFILE_2 or die "Cannot close filehandle to paired-end-2 file $opts{'output2'}: $!";

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

