#!/usr/bin/env perl

# WBpfm_2_WBpwm.pl -- Erich Schwarz <emsch@its.caltech.edu>, 6/12/2009; significantly revised, 10/10/2012.
# Purpose: convert .ace PFM to rough PWM; or, update PWM numbers from a PFM; allow rounding-off of PWM values.

use strict;
use warnings;
use Getopt::Long;
use TFBS::Matrix::PFM;  # From: http://tfbs.genereg.net/

my $text_source = q{};
my $num_source  = q{};
my $TEXT;
my $NUM;
my $clean;
my $help;

# Use this variable to allow rounding with user-specified precision:
my $DIGITS = 'no';

GetOptions ( 'text=s'          => \$text_source,
             'numbers=s'       => \$num_source,
             'clean'           => \$clean,
             'rounding=i{0,1}' => \$DIGITS, 
             'help'            => \$help
           );

# If user just says '-r ', then default to four digits.
$DIGITS or ($DIGITS = 4);
# If it's still 'no', *then* set to q{}:
if ($DIGITS eq 'no') { 
    $DIGITS = q{};
}

if ($help) { 
    die_loudly();
}

### Choose source(s) of data: ###

if ($text_source) { 
    if (! $num_source ) { 
        die "Failed to specify a source file for numbers!\n";
    }
    open $TEXT, '<', $text_source 
        or die "Can't open text source $text_source: $!";
}

if ($num_source) { 
    if (! $text_source ) {
        die "Failed to specify a source file for text!\n";
    }
    open $NUM, '<', $num_source 
        or die "Can't open number source file $num_source: $!";
}

if ( (! $text_source ) and (! $num_source ) and (! $TEXT ) and (! $NUM ) ) { 
    $TEXT = \*ARGV;
    $NUM  = \*ARGV;
}

### Uses of the script: ###

# Check to see if importing/printing leaves a PFM or PWM unchanged.
# Also a handy way to standardize file formats.
if ($clean) {
    my $pfms_ref = import_wbmat_linestream($TEXT);
    print_all_wbmats($pfms_ref);
    exit;
}

# Default usage, to map PFM file to raw PWM file:
if ( (! $clean) and ( $text_source eq $num_source ) ) { 
    my $pfms_ref = import_wbmat_linestream($TEXT);
    my $pwms_ref = pwmize_all_wb_pfms($pfms_ref,$DIGITS);
    print_all_wbmats($pwms_ref);
    exit;
}

# Alternative usage, to splice PFM-derived numbers into a PWM file:
if ( (! $clean ) and (-e $text_source) and (-e $num_source) and ( $text_source ne $num_source ) ) { 
    my $pwms_text_ref   = import_wbmat_linestream($TEXT);
    my $parent_pfms_ref = import_wbmat_linestream($NUM);
    my $raw_pwm_nos_ref = pwmize_all_wb_pfms($parent_pfms_ref,$DIGITS);
    my $new_pwms_ref    = transfer_pwm_nos($pwms_text_ref, $raw_pwm_nos_ref);
    print_all_wbmats($new_pwms_ref);
    exit;
}

### Subroutines: ### 

# Requires an incoming line stream.
sub import_wbmat_linestream { 
    my $_TEXT = $_[0];
    my $mats_ref;
    my $matrix = q{};

    while ( my $input = <$_TEXT> ) {
        chomp $input;

        # Get and remember identity of input matrix; reject duplex entries.
        if ( $input =~ / \A 
                         Position_Matrix 
                         \s+ : \s+ \" 
                         ( [^\"]+ ) 
                         \"  
                       /xms ) {

            $matrix = $1;
            if ( $mats_ref->{$matrix} ) {
                die "Error: matrix $matrix is being entered twice.\n";
            }
            $mats_ref->{$matrix}->{'Position_Matrix_line'} = $input;
        }

        # Require a single authorized matrix type without quotes.
        if ( $input =~ / \A
                          Type
                          \s+ 
                          ( \S+ )
                       /xms ) { 

            my $matrix_type = $1;
            if ( $matrix_type !~ /\A (Frequency|Weight) \z/xms ) {
                if ( $matrix_type =~ /(\"|\')/xms ) {
                    warn "Remove quotes from $matrix_type!\n";
                }
                die "Matrix $matrix is unauthorized type $matrix_type\n";
            }
            if ( $mats_ref->{$matrix}->{'type'} ) {
                die "Matrix $matrix has nonunique $matrix_type type\n";
            }
            $mats_ref->{$matrix}->{'type'}      = $matrix_type;
            $mats_ref->{$matrix}->{'Type_line'} = $input;
        }

        # Mult. Description/Remarks allowed, in ordered nonredund. array.
        my $allowed_Nliner = '(Brief_id|Consensus|Description|Remark)';
        if ( $input =~ / \A        
                         ($allowed_Nliner)
                         \s+ \" 
                         [^\"]+
                         \" 
                       /xms ) {

            my $Nline_tag = $1;

            if ( !$mats_ref->{$matrix}->{'seen'}->{$input} ) {
                push @{ $mats_ref->{$matrix}->{$Nline_tag} }, $input;
                $mats_ref->{$matrix}->{'seen'}->{$input} = 1;
            }
        }

        if ( $input =~ /  \A 
                          Site_values 
                          \s+ 
                          ([ACGT]) 
                          ( (\s+ \S+){1,} ) 
                          \s*
                       /xms ) {

            my $residue     = $1;
            my $number_line = $2;
            $number_line =~ s/\A\s+//xm;
            $mats_ref->{$matrix}->{'Site_values'}->{$residue} 
              = $number_line;
        }
        if ( $input =~ /  \A
                          Sites_used
                          \s+
                          \d+
                          \s*
                       /xms ) {
        
            my $sites_used = $1;
            $mats_ref->{$matrix}->{'Sites_used'} = $input;
        }   

        if ( $input =~ / \A 
                         (Derived_from_matrix
                         \s+
                         (\S+))
                         \s*
                       /xms ) {
            my $parent_txt    = $1;
            my $parent_matrix = $2;
            $mats_ref->{$matrix}->{'derived'}             = $parent_txt;
            $mats_ref->{$matrix}->{'parent_matrix'}       = $parent_matrix;
            # This can generate near-empty, non-print-worthy records:
            $mats_ref->{$parent_matrix}->{'child_matrix'} = $matrix;
        }
    }
    return $mats_ref;
}

# Input: $wb_pfms_ref.  Output: raw PWM file data.
sub pwmize_all_wb_pfms {
    my $pfms_ref = $_[0];
    my $_DIGITS  = $_[1];
    my $pwms_ref;

    foreach my $pfm ( sort keys %{$pfms_ref} ) {
        if ( $pfms_ref->{$pfm}->{'type'} ne 'Frequency' ) {
            die "Cannot convert non-frequency matrix $pfm to PWM.\n";
        }
        if ( $pfms_ref->{$pfm}->{'type'} eq 'Frequency' ) {
            my $pwm = $pfm . '.pwm';

            $pwms_ref->{$pwm}->{'type'}          = 'Weight';
            $pwms_ref->{$pwm}->{'parent_matrix'} = $pfm;
            $pwms_ref->{$pwm}->{'derived'}       = "Derived_from_matrix  $pfm";

            $pwms_ref->{$pwm}->{'Type_line'} 
                = $pfms_ref->{$pfm}->{'Type_line'};
            $pwms_ref->{$pwm}->{'Type_line'} =~ s/Frequency/Weight/xm;

            if ( exists $pfms_ref->{$pfm}->{'Description'} ) {
                @{ $pwms_ref->{$pwm}->{'Description'} } 
                  = @{ $pfms_ref->{$pfm}->{'Description'} };
            }

            if ( exists $pfms_ref->{$pfm}->{'Brief_id'} ) {
                @{ $pwms_ref->{$pwm}->{'Brief_id'} }
                  = @{ $pfms_ref->{$pfm}->{'Brief_id'} };
            }

            if ( exists $pfms_ref->{$pfm}->{'Consensus'} ) {
                @{ $pwms_ref->{$pwm}->{'Consensus'} }
                  = @{ $pfms_ref->{$pfm}->{'Consensus'} };
            }

            if ( exists $pfms_ref->{$pfm}->{'Sites_used'} ) {
                $pwms_ref->{$pwm}->{'Sites_used'} = $pfms_ref->{$pfm}->{'Sites_used'} ;
            }

            # Extract matrix-string counts from PFMs:
            my @residues 
                = sort keys %{ $pfms_ref->{$pfm}->{'Site_values'} };
            my $matrixstring = q{};

            # Implicitly sort by A, C, G, T:
            foreach my $residue (@residues) { 
                $matrixstring 
                    .= $pfms_ref->{$pfm}->{'Site_values'}->{$residue} 
                       . "\n";
            }

            # This next section must be correct, or the PWMs will be junk.
            # 
            # I'm hard-coding C. elegans nucleotide freqs.
            #     from WS200:
            # 
            #     a    32367419
            #     c    17780761
            #     g    17756942
            #     t    32367086
            # 
            # This works for C. elegans but becomes increasingly dubious as
            #     one applies this to more remote species.  If and when this 
            #     script is used for non-elegans sequences, then the script 
            #     should be rewritten to allow argument-provided residue frequencies.

            my $a_freq = ( 32_367_419 / 100_272_208 );
            my $c_freq = ( 17_780_761 / 100_272_208 );
            my $g_freq = ( 17_756_942 / 100_272_208 );
            my $t_freq = ( 32_367_086 / 100_272_208 );

            my $freqs_ref = { A => $a_freq,
                              C => $c_freq,
                              G => $g_freq,
                              T => $t_freq, };
            my %input_values = ( -matrixstring     => $matrixstring,
                                 -name             => $pwm,
                                 -ID               => $pwm,
                                 -bg_probabilities => $freqs_ref, 
                               );

            # Convert matrix strings from PFM to PWM values.
            my $tfbs_pfm  = TFBS::Matrix::PFM->new( %input_values ); 
            my $wb_tfbs_pwm = $tfbs_pfm->to_PWM();
            my $pwm_matrixstring = $wb_tfbs_pwm->rawprint();
            my %pwm_matrixstrings = ();
            @pwm_matrixstrings{ @residues } 
                = (split /\n/, $pwm_matrixstring);
            foreach my $residue (@residues) { 

                # TFBS gives absurdly precise decimals ...
                my $matrix_string = $pwm_matrixstrings{$residue};

                # So (as a conscious choice, not a default) round them off:
                if ($_DIGITS) { 
                    $matrix_string = round_off_pwm_no($matrix_string,$_DIGITS);
                }

                # Whether rounded or not:
                $pwms_ref->{$pwm}->{'Site_values'}->{$residue} 
                    = $matrix_string;
            }

            if ( exists $pfms_ref->{$pfm}->{'Remark'} ) {
                @{ $pwms_ref->{$pwm}->{'Remark'} } 
                    = @{ $pfms_ref->{$pfm}->{'Remark'} };
            }

            # N.B.: script fails if 'TFBS::Matrix::PFM->to_PWM()' is double-quoted:
            my $remark_header = 'Remark  "NOTE: weight matrix, derived by'
                              . ' TFBS::Matrix::PFM->to_PWM() from frequency matrix'
                              . " $pfm.\""
                              ;

            push @{ $pwms_ref->{$pwm}->{'Remark'} }, $remark_header;
        }
    }
    return $pwms_ref;
}

sub round_off_pwm_no { 
    my $input_matrixstring  = $_[0];
    my $output_matrixstring = q{};
    my $_DIGITS      = $_[1];

    while ( $input_matrixstring =~ /\A (\s*) (\S+) (.*) /xms ) { 
        my $leading_chars = $1;
        my $maybe_number  = $2;
        my $remainder     = $3;
        ( $maybe_number = sprintf("%.${_DIGITS}f", $maybe_number) ) 
            or die "Can't round off $maybe_number!\n";
        $output_matrixstring .= $leading_chars;
        $output_matrixstring .= $maybe_number;
        $input_matrixstring = $remainder;
    }
    # Append whatever's left:
    $output_matrixstring .= $input_matrixstring;
    return $output_matrixstring;
}

# Read and print data from $ref to *any* matrix (PFM or PWM) in .ace format.
sub print_all_wbmats {
    my $mats_ref = $_[0];
    print "\n";

    # Ignore near-empty matrix records, e.g., some parent-child links:
    my @printable = sort 
                    grep { exists($mats_ref->{$_}->{'Site_values'}) } 
                    keys %{ $mats_ref };

    foreach my $matrix ( @printable ) {

        print "Position_Matrix : \"$matrix\"\n";

        if ( exists $mats_ref->{$matrix}->{'Brief_id'} ) {
            my @brief_ids = @{ $mats_ref->{$matrix}->{'Brief_id'} };
            foreach my $brief_id (@brief_ids) {
                print "$brief_id\n";
            }
        }       

        if ( exists $mats_ref->{$matrix}->{'Description'} ) { 
            my @descrips = @{ $mats_ref->{$matrix}->{'Description'} };
            foreach my $descrip (@descrips) {
                print "$descrip\n";
            }
        }

        if ( exists $mats_ref->{$matrix}->{'Type_line'} ) {
            print $mats_ref->{$matrix}->{'Type_line'}, "\n";
        }

        if ( exists $mats_ref->{$matrix}->{'Consensus'} ) {
            my @consenses = @{ $mats_ref->{$matrix}->{'Consensus'} };
            foreach my $consensus (@consenses) {
                print "$consensus\n";
            }  
        }

        # This 'if' is kept as a failsafe:
        if ( exists $mats_ref->{$matrix}->{'Site_values'} ) {
            my @residues =
              sort keys %{ $mats_ref->{$matrix}->{'Site_values'} };
            foreach my $residue (@residues) {

                # Put this in to prevent additional spaces creeping in -- 
                #    a tiny thing which seems to be playing hob with 
                #    testing whether printing leaves data unaltered.
                my $vals_string = $mats_ref->{$matrix}->{'Site_values'}->{$residue};
                $vals_string =~ s/\A\s+//;
                $vals_string =~ s/\s+\z//;

                print 'Site_values  ', 
                      $residue, q{  },
                      $vals_string,
                      "\n",
                      ;
            }
        }

        if ( exists $mats_ref->{$matrix}->{'Sites_used'} ) {
            print "$mats_ref->{$matrix}->{'Sites_used'}\n";
        }

        if ( exists $mats_ref->{$matrix}->{'derived'} ) { 
            print "$mats_ref->{$matrix}->{'derived'}\n";
        }

        if ( exists $mats_ref->{$matrix}->{'Remark'} ) { 
            my @remarks = @{ $mats_ref->{$matrix}->{'Remark'} };
            foreach my $remark (@remarks) {
                print "$remark\n";
            }
        }
        print "\n";
    }
    return;
}

# Input: two refs., to PWMs with old text and new numbers.
# Output: single ref. to one PWM with both text and new numbers.
sub transfer_pwm_nos { 
    my ($pwms_ref, $pfm2pwm_nos_ref) = @_;
    foreach my $raw_pfm2pwm_matrix ( sort keys %{ $pfm2pwm_nos_ref } ) { 

        # For each PFM->PWM matrix (w/ nos.): see which PFM it came from.
        if (! $pfm2pwm_nos_ref->{$raw_pfm2pwm_matrix}->{'parent_matrix'} ) { 
            my @available_matrices = sort keys %{$pfm2pwm_nos_ref};
            warn "The only raw PFM->PWM matrices available are: @available_matrices\n";
            die "Raw PFM->PWM matrix $raw_pfm2pwm_matrix lacks well-identified PFM parent matrix!\n";
        }
        my $pfm_parent = $pfm2pwm_nos_ref->{$raw_pfm2pwm_matrix}->{'parent_matrix'};

        # For each PWM matrix (w/ text): use that PFM to elict its correct name.
        if (! $pwms_ref->{$pfm_parent}->{'child_matrix'} ) { 
            my @available_matrices = sort keys %{$pwms_ref};
            warn "The only PWM matrices available are: @available_matrices\n";
            die "PFM parent matrix $pfm_parent cannot be reliably mapped to child matrix!\n";
        }
        my $pwm_child = $pwms_ref->{$pfm_parent}->{'child_matrix'};

        my @residues = sort keys %{ $pfm2pwm_nos_ref->{$raw_pfm2pwm_matrix}->{'Site_values'} };
        foreach my $residue (@residues) { 
            if (! $pfm2pwm_nos_ref->{$raw_pfm2pwm_matrix}->{'Site_values'}->{$residue} ) { 
                die "No numerical data for residue $residue in raw PFM->PWM matrix $raw_pfm2pwm_matrix!\n";                
            }
            $pwms_ref->{$pwm_child}->{'Site_values'}->{$residue} 
                = $pfm2pwm_nos_ref->{$raw_pfm2pwm_matrix}->{'Site_values'}->{$residue};
        }
    }
    return $pwms_ref;
}

sub die_loudly {
    die "WBpfm_2_WBpwm.pl",
        " --text|-t [text source]",
        " --numbers|-n [number source]",
        " --clean|-c",
        " --rounding|-r [n, defaulting to 4 (digits)]",
        " --help|-h",
        "\n",
        ;
}

### Documentation (should perldoc): ###

# The default operation is:
#     Read in a PFM file.
#     Map it to a new, raw PWM file (which should not need hand-editing).
#     Print the PWM file.
#
# Note that the default works either if -t and -n name the same file,
#     or if they name nothing, and WBpfm_2_WBpwm.pl is used as
#     a naive UNIX-style filter!  (Since q{} eq q{}...)

# An alternative use of this script is:
#    Read in a PWM file; store its data (mainly, its non-numerical data).
#    Read in a file of PFM parents.
#    Generate new PWM numbers from the parents.
#    Copy these numbers to the first PWM files stored data, overwriting older numbers.
#    Print out the resulting PWM.
#       
# This was done to allow fast reliable corrections if (as happened way too many
#    times in 2008-2009) I stupidly got the PFM->PWM numerical mapping wrong.

# Example of older input, circa 2010:
#
# Position_Matrix : "WBPmat00000002"  // SKN-1.pfm
# Description    "SKN-1 homomeric binding site, frequency matrix." Paper_evidence "WBPaper00002041"
# Type           Frequency
# Site_values    A  6.25   17  14  7.25   14  0   0   33  0   1   13  12
# Site_values    C  12.25  4   2   3.25   0   0   33  0   0   15  12  7
# Site_values    G  3.25   2   3   1.25   19  0   0   0   0   9   2   2
# Site_values    T  11.25  10  14  21.25  0   33  0   0   33  8   6   12
# Remark         "Inferred by curator from mixed frequency data given in paper; threshold not clear from paper."

# Example of newer input, Oct. 2012:
#
