#!/usr/bin/env perl

# revcomp_fastq.pl -- Erich Schwarz <emsch@its.caltech.edu>, 1/29/2011.
# Purpose: properly reverse-complement seqs. from FASTQ files (and just reverse their quality lines); accept piped data via '-'.

use strict;
use warnings;
use Getopt::Long;

my $head1  = q{};
my $nt_seq = q{};
my $head2  = q{};
my $quals  = q{};

my %opts   = ();

GetOptions ( 'fastq=s' => \$opts{'input_fastq'}, 
             'help'    => \$opts{'help'},         );

if ( $opts{'help'} or (! $opts{'input_fastq'}) ) { 
    die "\n",
        "Format: revcomp_fastq.pl\n",
        "          --fastq|-f    [large FASTQ to reverse-complement (opt. '-')]\n",
        "          --help\n",
        "\n",
        ;
}

# Accept either a stream from '-' or a standard file.
my $INPUT_FASTQ;
if ($opts{'input_fastq'} eq '-') {
    # Special case: get the stdin handle
    $INPUT_FASTQ = *STDIN{IO};
}
else {
    # Standard case: open the file
    open $INPUT_FASTQ, '<', $opts{'input_fastq'} or die "Can't open $opts{'input_fastq'}. $!\n";
}

while (<$INPUT_FASTQ>) {
    $head1 = $_;
    $nt_seq = <$INPUT_FASTQ>;
    $head2  = <$INPUT_FASTQ>;
    $quals  = <$INPUT_FASTQ>;

    chomp( $head1, $nt_seq, $head2, $quals );
    if ( ( $head1 !~ /\A @ \S+ /xms ) or ( $head2 !~ /\A \+ /xms ) ) { 
        warn "Can't parse one or both of these headers:\n";
        warn "$head1\n";
        warn "$head2\n";
        die;
    }

    $nt_seq = flipseq( $nt_seq, 'revcomp' );
    $quals  = flipseq( $quals,  'reverse' );    
    print "$head1\n", "$nt_seq\n", "$head2\n", "$quals\n", ;
}
close $INPUT_FASTQ or die "Can't close filehandle to FASTA file $opts{'input_fastq'}: $!";

# Note that -- appropriately -- this does nothing to 'n' or 'N' residues.
sub flipseq { 
    my $in_string = $_[0];
    my $option    = $_[1];
    my %OK_opts = ( 'rev'     => 1,
                    'reverse' => 1,
                    'revcomp' => 1, );
    if (! exists $OK_opts{$option} ) {
        die "sub flipseq can't parse option \"$option\" for input text \"$in_string\"!\n";
    }

    if ( $option eq 'revcomp' ) { 
        $in_string =~ tr/[acgtACGT]/[tgcaTGCA]/;
    }
    my @in_residues = split //, $in_string;
    @in_residues = reverse @in_residues;
    my $out_string = join q{}, @in_residues;
    return $out_string;
}

