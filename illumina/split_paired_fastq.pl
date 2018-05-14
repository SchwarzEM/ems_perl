#!/usr/bin/env perl

# split_paired_fastq.pl -- Erich Schwarz <emsch@its.caltech.edu>, 1/27/2011.
# Purpose:  Take four lines at a time.  Check for superficially correct headers.  Split into two end-1 and end-2 files.

use strict;
use warnings;
use Getopt::Long;

my $input = q{};

my $head1  = q{};
my $nt_seq = q{};
my $head2  = q{};
my $quals  = q{};

my %opts  = ();

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
    "help"    => \$opts{'help'},
);

if ( (! $opts{'input'}) or (! $opts{'r1'}) or (! $opts{'r2'}) ) {
    $opts{'help'} = 1;
}

if ( $opts{'help'} ) {
    print "usage: split_paired_fastq.pl\n";
    print "       -i|--input  <name1>    input file name (fastq), mandatory (though '-' will work)\n";
    print "       --o1        <name2>    output paired-end 1 file name (fastq); default \"name1.output1\"\n";
    print "       --o2        <name3>    output paired-end 2 file name (fastq); default \"name1.output2\"\n";
    print "       --r1        <suffix1>  suffix marking paired-end read 1, default \"#0\\/1\".\n";
    print "       --r2        <suffix2>  suffix marking paired-end read 2, default \"#0\\/2\".\n";
    print "       -h|--help              help - print this message\n";
    exit;
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

open my $OUTFILE_1, '>', $opts{'output1'} or die "Cannot open output paired-end-1 file $opts{'output1'}: $!";
open my $OUTFILE_2, '>', $opts{'output2'} or die "Cannot open output paired-end-2 file $opts{'output2'}: $!";

while (<$INFILE>) { 
    # Tricky Perl syntax: 
    #     it forces me to start with $_ already loaded from <$INFILE>,
    #     but afterwards needs prompting for linefeeds.
    $head1 = $_;
    chomp $head1;
    $nt_seq = <$INFILE>;
    chomp $nt_seq;
    $head2  = <$INFILE>;
    chomp $head2;
    $quals  = <$INFILE>;
    chomp $quals;

    # Note: for reasons that are highly unclear to me as of 10/18/2010, the normally-good
    #    'xms' modifier on regular expressions causes parsing of header names from their
    #    paired-end suffixes to totally fail!
    #    So, to get the script to simply work, primitive non-legible non-xms regexes are 
    #        used whenever FASTQ headers and suffixes need to be sorted out.

    if ( (       ( $head1 !~ /\A[@]\S+$opts{'r1'}\z/ ) 
             and ( $head1 !~ /\A[@]\S+$opts{'r2'}\z/ ) ) 
         or ( $head2 !~ /\A \+ /xms                    ) ) { 
        warn "Can't parse one or both of these headers:\n";
        warn "$head1\n";
        warn "$head2\n";
        die;
    }

    elsif (     ( $head1 !~ /\A[@]\S+$opts{'r1'}\z/ ) 
            and ( $head1 !~ /\A[@]\S+$opts{'r2'}\z/ ) ) { 
        die "Header $head1 somehow has both suffixes $opts{'r1'} and $opts{'r2'}!\n";
    }

    else { 
        if ( $head1 =~ /\A[@](\S+)$opts{'r1'}\z/ ) {
            my $print_text = join "\n", ( $head1, $nt_seq, $head2, $quals, );
            print $OUTFILE_1 "$print_text\n";
        }
        elsif ( $head1 =~ /\A[@](\S+)$opts{'r2'}\z/ ) {
            my $print_text = join "\n", ( $head1, $nt_seq, $head2, $quals, ); 
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

