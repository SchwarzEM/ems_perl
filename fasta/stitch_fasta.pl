#!/usr/bin/env perl

# stitch_fasta.pl -- Erich Schwarz <emsch@its.caltech.edu>, originally 3/14/2009; updated 7/11/2016.
# Purpose: given user-specified orders, convert FASTA to glued FASTA.

use strict;
use warnings;
use autodie;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

# We use line-commands to specify how many 'N' residues are used as glue; the default is 1 residue.
my $GLUE_SEQUENCE = 'N';

my $fasta = q{};
my $order = q{};
my $glue  = 1;
my $rest;

my @listed_seqs = ();
my $start_seq   = q{};
my %order_line  = ();

my $seq_name     = q{};
my %seqs2headers = ();
my %sequences    = ();

my $help;

GetOptions ( 'fasta=s' => \$fasta,
             'order=s' => \$order, 
             'glue:i'  => \$glue,
             'rest'    => \$rest,
             'help'    => \$help);

if ( (! -r $fasta ) or (! -r $order) or (! looks_like_number($glue) ) or $help) { 
    die "Format: stitch_fasta.pl\n",
        "    -f|--fasta [FASTA input]\n",
        "    -o|--order [file listing splice order input of FASTA members, one per line, as '+seq1', '-seq2', etc.]\n",
        "    -g|--glue  [number of 'N' residues to use as glue; default value is 1; must be positive integer\n",
        "    -r|--rest  [also print out any *un*spliced FASTA sequences as part of output; optional, not default\n",
        "    -h|--help  [print this message]\n",
        ;
}

if ( ( $glue != int($glue) ) or ( $glue <= 0 ) ) {
    die "--glue must be positive integer, not \"$glue\"\n";
}
$GLUE_SEQUENCE = ($GLUE_SEQUENCE x $glue);

# Read in the marching orders for what sequences to stitch and how.
open my $ORDER, '<', $order;
while (my $input = <$ORDER>) { 
    chomp $input;

    if ( $input =~ / \A [\+\-]\S+ \z /xms ) {
        push @listed_seqs, $input;
    }
    else {
        die "Couldn't parse input line \"$input\" from order file: $order!\n";
    }
}
close $ORDER;

# We will key the hash of jobs with the plain name of the first seq.:
$start_seq = $listed_seqs[0];
$start_seq =~ s/\A[\+|\-]//;

# Key by first name; store our marching orders as an array reference:
$order_line{$start_seq} = \@listed_seqs;

# Read the entire FASTA into a hash-table in RAM.  Works OK unless input data are huge.
open my $FASTA, '<', $fasta;
while (my $input = <$FASTA>) { 
    chomp $input;
    if ($input =~ /\A > ( (\S+) .*) /xms) { 
        $seq_name = $2; 
        $seqs2headers{$seq_name} = $1;
        $sequences{$seq_name} = q{};
    }
    elsif ($input =~ /[a-zA-Z]/) { 
        $input =~ s/[^a-zA-Z]//g;
        $sequences{$seq_name} .= $input;
    }
}
close $FASTA;

# Now, do a special retrieval and output for things with marching orders:
foreach my $seq_name2 (sort keys %order_line) { 
    # Defensive programming, to prevent lossage from bugs.
    reality_check($seq_name2);

    # Start with these empty, for each list of seqs. to stitch.
    my @aggregate_names   = ();
    my @aggregate_headers = ();
    my @output_seqs       = ();

    # Deference the job list which we stored as an array reference:
    my @listed_job_names = @{ $order_line{$seq_name2} } ;

    # OK, plonk through each stiching.
    foreach my $job (@listed_job_names) { 

        # Initialize:
        my $orientation = q{};
        my $next_seq = q{};

        # At last!
        if ($job =~ /\A ([\+|\-])(\S+) /xms ) {
            $orientation = $1;
            $next_seq = $2;
            reality_check($next_seq);
            if ($orientation eq '+') { 
                push @output_seqs, $sequences{$next_seq};
            }
            if ($orientation eq '-') {
                push @output_seqs, revcomp($sequences{$next_seq}); 
            }
            push @aggregate_names, $next_seq;
            push @aggregate_headers, "[$job -- $seqs2headers{$next_seq}]";
        }
        # Buggy inputs die here:
        else {
            die "Can't parse job $job!\n";
        }

        # Next step is CRUCIAL.  We destroy each record in the general FASTA hash data
        #    as we go.  This defeats any effort to read a sequence twice.
        #    It also makes the later general read-out loop (optionally invoked with 'rest') safe for $seq_name3, 
        #    since only normal unspliced sequences have any chance of being read.

        delete $sequences{$next_seq};
        delete $seqs2headers{$next_seq};
    }

    # Now get ready to spit out a FASTA record of a spliced sequence set:
    my $final_name     = join '_', @aggregate_names;
    my $final_header   = join '; ', @aggregate_headers;
    my $final_sequence = join "$GLUE_SEQUENCE", @output_seqs;

    print '>', $final_name, ' ', $final_header, "\n";
    regular_FASTA_seq_print($final_sequence);
}

# Optional: also print out the whole remaining list of any *unspliced* sequences.
if ($rest) {
    foreach my $seq_name3 (sort keys %sequences) { 
        print '>', "$seqs2headers{$seq_name3}\n";
        regular_FASTA_seq_print($sequences{$seq_name3});
    }
}

sub revcomp { 
    my $input_sequence = $_[0];
    $input_sequence = reverse($input_sequence);
    if ( $input_sequence =~ /[^acgtnACGTN]/xms ) { 
        die "Not currently designed to parse any",
            " letters but a, c, g, t, n, A, C, G, T, or N!\n",
            ;
    }
    $input_sequence =~ tr/acgtnACGTN/tgcanTGCAN/;
    return $input_sequence;
}

sub reality_check { 
    my $_query_seq = $_[0];
    if ( (! $sequences{$_query_seq} ) or (! $seqs2headers{$_query_seq} ) ) {
        die "Failed to find $_query_seq in input FASTA data!\n";
    }
    # Pro forma:
    return;
}

sub regular_FASTA_seq_print { 
    my $raw_seq_text = $_[0];
    my @output_lines
        = unpack("a60" x (length($raw_seq_text)/60 + 1), $raw_seq_text);
    foreach my $output_line (@output_lines) {
        if ($output_line =~ /\S/) {
            print "$output_line\n";
        }
    }   
    # Pro forma:
    return;
}

