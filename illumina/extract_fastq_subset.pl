#!/usr/bin/env perl

# extract_fastq_subset.pl -- Erich Schwarz <emsch@its.caltech.edu>, 11/2/2010.
# Purpose: given seqlist file or single-name argument, either extract or exclude from FASTQ files; accept piped data via '-'; warn of misses.

use strict;
use warnings;
use Getopt::Long;

my $head1  = q{};
my $nt_seq = q{};
my $head2  = q{};
my $quals  = q{};

my %opts   = ();

GetOptions ( 'list=s'  => \$opts{'input_list'},
             'fastq=s' => \$opts{'input_fastq'}, 
             'exclude' => \$opts{'exclude'},
             'help'    => \$opts{'help'},         );

if ( $opts{'help'} or (! $opts{'input_list'}) or (! $opts{'input_fastq'}) ) { 
    die "\n",
        "Format: extract_fastq_subset.pl\n",
        "          --list|-l     [seqs. list OR seq. name]\n",
        "          --fastq|-f    [large FASTQ to extract (opt. '-')]\n",
        "          --exclude|-e  [optional -- reject list, keep rest]\n",
        "          --help\n",
        "\n",
        ;
}

my %input_names    = ();

my $warnings      = $opts{'input_fastq'} 
                    . q{.} 
                    . $opts{'input_fastq'} 
                    . '.warnings'
                    ;

if (-e $opts{'input_list'}) { 
    open my $INPUT_LIST, "$opts{'input_list'}" or die "Can't open input namelist $opts{'input_list'}. $!\n";
    while (my $inline = <$INPUT_LIST>) { 
        chomp $inline;
        if ( $inline =~ /\A \s* (\S+) \s* /xms) { 
            my $seq = $1;
            $input_names{$seq} = 1;
        }
    }
    close $INPUT_LIST;
}

# If no list file, treat as single-name argument.
if (!-e $opts{'input_list'}) { 
    $input_names{ $opts{'input_list'} } = 1;
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

    if ( ( $head1 !~ /\A @ \S+ /xms ) or ( $head2 !~ /\A \+ /xms ) ) { 
        warn "Can't parse one or both of these headers:\n";
        warn "$head1\n";
        warn "$head2\n";
        die;
    }

    if ( ( $head1 =~ / \A @ ( \S+ ) \s* /xms ) 
         and (    ( $input_names{$1} and (! $opts{'exclude'})          ) 
               or ( $opts{'exclude'} and (! $input_names{$1} ) ) 
             ) 
       ) { 
        my $seq = $1;
        print $head1, $nt_seq, $head2, $quals, ;
        delete $input_names{$seq};  # not 'undef'; Perl Cookbook 5.4.
    }
}
close $INPUT_FASTQ or die "Can't close filehandle to FASTA file $opts{'input_fastq'}: $!";

# If this is a 'normal' extraction, not an exclusion-filtering:
if (! $opts{'exclude'} ) { 
    # Scan hash for any non-zero values:
    my $missing_count = scalar (keys %input_names);

    # If any, warn they were missed:
    if ( $missing_count >= 1) { 
        open my $WARNINGS, '>', $warnings 
            or die "Can't open warnings file $warnings. $!";
        foreach my $key (sort keys %input_names) {
            if ( exists $input_names{$key} ) {   # not 'defined'; Perl Cook. 5.2.
                print {$WARNINGS} "$key not found in $opts{'input_fastq'}\n";
            }
        }
        close $WARNINGS 
            or die "Can't close filehandle to warnings file $warnings: $!";
    }
}

