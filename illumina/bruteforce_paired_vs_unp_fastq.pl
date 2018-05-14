#!/usr/bin/env perl

# bruteforce_paired_vs_unp_fastq.pl -- Erich Schwarz <emsch@its.caltech.edu>, 10/18/2010.
# Purpose:  Take four lines at a time.  Check for superficially correct headers.  Fully namesort (with much RAM and brute force!) into truly paired+sorted and unpaired output files.

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

# Default suffixes for paired-end read 1 and paired-end read 2:
$opts{'r1'} = '#0\/1';
$opts{'r2'} = '#0\/2';

GetOptions(
    "input=s"           => \$opts{'input'},
    "paired_output:s"   => \$opts{'paired'},  
    "unpaired_output:s" => \$opts{'unpaired'},
    "r1=s"              => \$opts{'r1'},
    "r2=i"              => \$opts{'r2'}, 
    "help"              => \$opts{'help'},
);

if ( (! $opts{'input'}) or (! $opts{'r1'}) or (! $opts{'r2'}) ) {
    $opts{'help'} = 1;
}

if ( $opts{'help'} ) {
    print "usage: paired_vs_unp_fastq.pl\n";
    print "       -i|--input            <name1>    input file name (fastq), mandatory\n";
    print "       -p|--paired_output    <name2>    output paired-end file name (fastq); default \"name1.paired\"\n";
    print "       -u|--unpaired_output  <name3>    output unpaired-end file name (fastq); dault \"name1.unpaired\"\n";
    print "       --r1                  <suffix1>  suffix marking paired-end read 1, default \"#0/1\".\n";
    print "       --r2                  <suffix1>  suffix marking paired-end read 2, default \"#0/1\".\n";
    print "       -h|--help                        help - print this message\n";
    exit;
}

if (! $opts{'paired'} ) { 
    $opts{'paired'} = $opts{'input'} . ".paired";
    $opts{'paired'} = safename( $opts{'paired'} );
}

if (! $opts{'unpaired'} ) {
    $opts{'unpaired'} = $opts{'input'} . ".unpaired";
    $opts{'unpaired'} = safename( $opts{'unpaired'} );
}

open my $INFILE,   '<', $opts{'input'}    or die "Cannot open input jumbled-reads file $opts{'input'}: $!";
open my $PAIRED,   '>', $opts{'paired'}   or die "Cannot open output paired-reads file $opts{'paired'}: $!";
open my $UNPAIRED, '>', $opts{'unpaired'} or die "Cannot open output unpaired-red file $opts{'unpaired'}: $!";

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
        # Overspecify this to avoid confusing the Perl interpreter.
        my @read_data     = ( $head1, $nt_seq, $head2, $quals, );
        my $read_data_ref = \@read_data;

        if ( $head1 =~ /\A[@](\S+)$opts{'r1'}\z/ ) {
            $read_stem = $1;
            if ( exists $data_ref->{$read_stem}->{'r1'} ) {
                warn "Multiple sequence entries with header $head1; overwriting data!\n";
            }
            $data_ref->{$read_stem}->{'r1'} = $read_data_ref;
        }
        elsif ( $head1 =~ /\A[@](\S+)$opts{'r2'}\z/ ) {
            $read_stem = $1;
            if ( exists $data_ref->{$read_stem}->{'r2'} ) { 
                warn "Multiple sequence entries with header $head1; overwriting data!\n";
            }
            $data_ref->{$read_stem}->{'r2'} = $read_data_ref;
        } 
    }
}
close $INFILE, or die "Cannot open close filehandle to jumbled-reads file $opts{'input'}: $!";

my @read_stems = sort keys %{ $data_ref };
foreach my $read_stem2 (@read_stems) { 
    # Set to print to the paired file, if both ends of a read are available:
    if ( ( exists $data_ref->{$read_stem2}->{'r1'} ) and ( exists $data_ref->{$read_stem2}->{'r2'} ) ) {
        select $PAIRED;
    }
    # Set to print to the unpaired file, otherwise:
    if (    ( ( exists $data_ref->{$read_stem2}->{'r1'} ) and (! exists $data_ref->{$read_stem2}->{'r2'} ) ) 
         or ( (! exists $data_ref->{$read_stem2}->{'r1'} ) and ( exists $data_ref->{$read_stem2}->{'r2'} ) ) ) {
        select $UNPAIRED;
    }
    # Print whatever exists to whatever's been selected:
    if ( exists $data_ref->{$read_stem2}->{'r1'} ) { 
        my $print_text = join "\n", @{ $data_ref->{$read_stem2}->{'r1'} };
        print "$print_text\n";
    }
    if ( exists $data_ref->{$read_stem2}->{'r2'} ) { 
        my $print_text = join "\n", @{ $data_ref->{$read_stem2}->{'r2'} };
        print "$print_text\n";
    }
    # Restore print etc. output to default!
    select STDOUT; 
}

close $PAIRED   or die "Cannot close filehandle to paired-reads file $opts{'paired'}: $!";
close $UNPAIRED or die "Cannot close filehandle to paired-reads file $opts{'paired'}: $!";

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

