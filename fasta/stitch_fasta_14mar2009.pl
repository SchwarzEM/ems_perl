#!/usr/bin/env perl

# stitch_fasta_14mar2009.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/14/2009.
# Purpose: given user-specified orders, convert FASTA to glued FASTA.  [Legacy version; kept in case I need to rerun older work.]

use strict;
use warnings;
use Getopt::Long;

# Hard-coded to stitch contigs with just one 'N'; could be made more flexible!
my $GLUE_SEQUENCE = 'N';

my $fasta = q{};
my $order = q{};
my %order_line = ();

my $seq_name     = q{};
my %seqs2headers = ();
my %sequences    = ();

GetOptions ( 'fasta=s' => \$fasta,
             'order=s' => \$order, );

if ( (! -r $fasta ) or (! -r $order) ) { 
    die "Format: ./stitch_fasta.pl",
        " -f|--fasta [FASTA input]",
        " -o|--order [splice order input]\n",
        ;
}

# Read in the marching orders for what sequences to stitch and how.
open my $ORDER, '<', $order 
    or die "Can't open stitch-order file $order: $!";
while (my $input = <$ORDER>) { 
    chomp $input;

    # Make life easy later -- die loudly if bad syntax.
    if ( $input !~ / \A [\+\-]\S+ ( , \s+ [\+|\-]\S+ )+ \s* \z /xms ) { 
        die "Couldn't parse input line \"$input\" from order file: $order!\n";
    }

    # If good syntax, record the orders for later.
    if ( $input =~ / \A ( [\+\-]\S+ ( , \s+ [\+|\-]\S+ )+ ) \s* \z /xms ) {
        $input = $1;

        # Tolerate, but delete, trailing commas.
        $input =~ s/\,\z//;

        my @listed_seqs = split /,\s+/, $input;

        # We will key the hash of jobs with the plain name of the first seq.:
        my $start_seq = $listed_seqs[0];
        $start_seq =~ s/\A[\+|\-]//;

        # Key by first name; store our marching orders as an array reference:
        
        $order_line{$start_seq} = \@listed_seqs;
    }
}
close $ORDER or die "Can't close filehandle to stitch-order file $order: $!";

# Read the entire FASTA into a hash-table in RAM.  Works OK unless data huge...
open my $FASTA, '<', $fasta 
    or die "Can't open FASTA file $fasta: $!";

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
close $FASTA or die "Can't close filehandle to FASTA file $fasta: $!";

# Now, do a special retrieval and output for things with marching orders:
foreach my $seq_name2 (sort keys %order_line) { 

    # Very 'defensive' programming, but should prevent lossage from bugs.
    reality_check($seq_name2);

    # Start with these empty, for each list of seqs. to stitch.
    my @aggregate_names   = ();
    my @aggregate_headers = ();
    my @output_seqs = ();

    # Deference the job list which we stored as an array reference:
    my @listed_job_names = @{ $order_line{$seq_name2} } ;

    # OK, plonk through each stiching.
    foreach my $job (@listed_job_names) { 

        # Initialize:
        my $orientation = q{};
        my $next_seq = q{};

        # Again, very defensive, but buggy inputs die here:
        if ($job !~ /\A [\+|\-]\S+ /xms ) { 
            die "Can't parse job $job!\n";
        }

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

        # Next step is CRUCIAL.  We destroy each record in the general FASTA hash data
        #    as we go.  This defeats any effort to read a sequence twice.
        #    It also makes the later general read-out loop safe for $seq_name3, 
        #    since only normal unspliced sequences have any chance of being read.

        delete $sequences{$next_seq};
        delete $seqs2headers{$next_seq};
    }

    # Now get ready to spit out a FASTA record of a spliced sequence set:
    my $final_name = join '_', @aggregate_names;
    my $final_header = join '; ', @aggregate_headers;
    my $final_sequence = join "$GLUE_SEQUENCE", @output_seqs;

    print '>', $final_name, ' ', $final_header, "\n";
    regular_FASTA_seq_print($final_sequence);
}

foreach my $seq_name3 (sort keys %sequences) { 
    print '>', "$seqs2headers{$seq_name3}\n";
    regular_FASTA_seq_print($sequences{$seq_name3});
}

sub revcomp { 
    my $input_sequence = $_[0];
    $input_sequence = reverse($input_sequence);
    if ( $input_sequence =~ /[^acgtACGT]/xms ) { 
        die "Not currently designed to parse any",
            " letters but a, c, g, t, A, C, G, or T!\n",
            ;
    }
    $input_sequence =~ tr/acgtACGT/tgcaTGCA/;
    return $input_sequence;
}

sub reality_check { 
    my $query_seq = $_[0];
    if ( (! $sequences{$query_seq} ) or (! $seqs2headers{$query_seq} ) ) {
        die "Failed to find $query_seq in input FASTA data!\n";
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

