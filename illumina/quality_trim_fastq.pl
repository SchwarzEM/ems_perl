#!/usr/bin/env perl

# quality_trim_fastq.pl -- Erich Schwarz <ems394@cornell.edu>, 11/19/2019.
# Purpose:  Take four lines at a time.  Print only, as chosen: all but last (unreliable) residue; residues up to maximum length; residues up to last non-N; residues >= quality threshold; residues >= minimum length.  NOTE: no longer *requires* the user to specify the offset X for Phred+X -- assumes default 33 ; X was 64 for Illumina reads until <= 10/2011 but then became 33; at some point in the more distant past, it was something else (neither 64 nor 33...)

use strict;
use warnings;
use Getopt::Long;

my $input  = q{};

my $head1  = q{};
my $nt_seq = q{};
my $head2  = q{};
my $quals  = q{};

my %opts   = ();

# Default is to allow any amount of trimming, but this can lead to vacant sequence lines!
$opts{'minimum'} = 0;

# Default quality offset for FastQ is 33; can be overridden by user, for dealing with old / weird data.
$opts{'quality'} = 33;

GetOptions(
    "input=s"     => \$opts{'input'},
    "output=s"    => \$opts{'output'},
    "quality=i"   => \$opts{'quality'},

    "end_trim:i"  => \$opts{'end_trim'},
    "uppermost:i" => \$opts{'maximum'},
    "no_Ns"       => \$opts{'no_Ns'},
    "threshold:i" => \$opts{'threshold'}, 
    "minimum:i"   => \$opts{'minimum'}, 

    "chastity"    => \$opts{'chastity'},

    "help"        => \$opts{'help'},
);

if (    (! $opts{'input'}             ) 
     or (! $opts{'output'}            )
     or (! $opts{'quality'}           ) 
     or (! (    $opts{'end_trim'} 
             or $opts{'no_Ns'} 
             or $opts{'threshold'} 
             or $opts{'minimum'} 
             or $opts{'maximum'}    ) ) ) {
    $opts{'help'} = 1;
}

if ( $opts{'help'} ) {
    print "\n";
    print "usage: quality_trim_fastq.pl\n";
    print "\n";
    print "       -i|--input      <in>       input file name (fastq); opt. stream input via '-'\n";
    print "       -o|--output     <out>      output file name (fastq), required\n";
    print "\n";
    print "       [at least one of the following arguments is required]\n";
    print "       -e|--end_trim   <integer>  1. Before anything else, trim this many residues off end (typically, 1 last unreliable residue)\n";
    print "       -u|--uppermost  <integer>  2. Maximum residue length - trim all reads to this length, immediately after --end_trim\n";
    print "       -n|--no_Ns                 3. Trim all *starting* N residues, and then all 3' residues to last non-N\n";
    print "       -t|--threshold  <integer>  4. Quality threshold - trim all 3'-ward residues below this score\n";
    print "       -m|--minimum    <integer>  5. Minimum residue length - after all other filters, censor any reads trimmed below this length\n";
    print "\n";
    print "       [the following value, if not specified by user, is 33 (Phred+33) by default]\n";
    print "       -q|--quality    <integer>  quality score offset for FASTQ: for Illumina, 33 (Phred+33) in Oct. 2011+; 64 (Phred+64) in 2009-2011\n";
    print "\n";
    print "       [the following option is allowed, but not recommended]\n";
    print "       -c|--chastity              6. Override the automatic filtering out of any read with ", '1:Y:\d or 2:Y:\d', " in the header\n", ;
    print "\n";
    print "       [print this message:]\n";
    print "       -h|--help\n";
    print "\n";
    exit;
}


# Accept either a stream from '-' or a standard file.
my $INFILE;
if ($opts{'input'} eq '-') {
    # Special case: get the stdin handle
    $INFILE = *STDIN{IO};
}
else {
    # Standard case: open the file
    open $INFILE,  '<', $opts{'input'} or die "Cannot open input file $opts{'input'}: $!";
}

open my $OUTFILE, '>', $opts{'output'} or die "Cannot open output file $opts{'output'}: $!";

# If no positive-integer minimum length given, make it zero:
if ( (! $opts{'minimum'} ) or ( $opts{'minimum'} < 1 ) ) {
    if (! ( $opts{'minimum'} == 0 ) ) { 
        warn "Not able to accept minimum length of $opts{'minimum'}; resetting to 0!\n";
        $opts{'minimum'} = 0;
    }
}

# if no positive-integer maximum length given, delete entirely:
if ( $opts{'maximum'} and ( $opts{'maximum'} < 1 ) ) { 
    warn "Not able to accept maximum length of $opts{'maximum'}; deleting entirely!\n";
    delete $opts{'maximum'};
}

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

    if ( ( $head1 !~ /\A @ \S+ /xms ) or ( $head2 !~ /\A \+ /xms ) ) { 
        warn "Can't parse one or both of these headers:\n";
        warn "$head1\n";
        warn "$head2\n";
        die;
    }

    # Begin by allowing full length of existing read.
    my $final_length = length($nt_seq);

    # Trim it back by $opts{'end_trim'} residues:
    if ( $opts{'end_trim'} and ( $opts{'end_trim'} >= 1 ) ) { 
        $final_length = ($final_length - $opts{'end_trim'});
    }

    # Trim it back to $opts{'maximum'}, if that's been defined:
    if ( $opts{'maximum'} and ( $opts{'maximum'} < $final_length ) ) { 
        $final_length = $opts{'maximum'};
    }

    # Trim it back to last non-N residue:
    if ( $opts{'no_Ns'} and ( $nt_seq =~ /[nN]/xms ) ) {
        my $non_N_nt_seq = $nt_seq;

        # This is crucial to get right, or Bad Things will happen:
        #     lines will be completely deleted, not just freed of leading or 3'-ward N residues, 
        #     or single embedded N residues will be left intact!

        my $lead_Ns_seq = q{};
        my $lead_Ns_len = 0;
        if ( $non_N_nt_seq =~ /\A ([nN]+) .* \z /xms ) { 
            # Count leading Ns:
            $lead_Ns_seq  = $1;
            $lead_Ns_len  = length($lead_Ns_seq);

            # Use precise substr to enforce correct 5'-ward truncation of both seq. and quals.
            $nt_seq       = substr($nt_seq, $lead_Ns_len);
            $non_N_nt_seq = substr($non_N_nt_seq, $lead_Ns_len);
            $quals        = substr($quals,  $lead_Ns_len);
        }
        
        # Then delete all residues from the first non-leading N residue, onwards:
        $non_N_nt_seq =~ s/[nN].*\z//;

        # Adamantly enforce successful editing.
        if ( ( $non_N_nt_seq !~ /\A [ACGTacgt]* \z/xms ) or ( $non_N_nt_seq =~ /[nN]/xms ) ) { 
            die "Can't correctly remove N (or Ns) from sequence (in read $head1): $nt_seq\n";
        }

        my $non_N_length = length($non_N_nt_seq);
        if ( $non_N_length < $final_length ) { 
            $final_length = $non_N_length;
        }
    }

    # Trim it back to last residue of a given quality:
    if ($opts{'threshold'}) {
        my $good_nt     = 0;
        my @qual_scores = split //, $quals;
        TRIM: foreach my $score (@qual_scores) { 
            $score = ord($score) - $opts{'quality'};  # Was once just '- 64'; but Illumina shifted back from Phred+64 to Phred+33!
            if ( $score >= $opts{'threshold'} ) { 
                $good_nt++;
            }
            else { 
                last TRIM;  # Stop and print trimmed reads.        
            }
        }
        if ( $good_nt < $final_length ) { 
            $final_length = $good_nt;
        }
    }

    # No point in doing anything unless $final_length >= $opts{'minimum'}.
    if (     ( $final_length >= $opts{'minimum'}                               ) 
         and ( $opts{'chastity'} or ( $head1 !~ /\A @ .* (?:1|2):Y:\d /xms )       ) ) {
        # Compute the desired sequence and quality strings:
        $nt_seq = substr($nt_seq, 0, $final_length);
        $quals  = substr($quals,  0, $final_length);

        # Print all four lines of the FASTQ entry:
        print $OUTFILE "$head1\n";
        print $OUTFILE "$nt_seq\n";
        print $OUTFILE "$head2\n";
        print $OUTFILE "$quals\n";
    }
}

