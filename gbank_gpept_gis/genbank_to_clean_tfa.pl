#!/usr/bin/perl -w

# Program: genbank_to_clean_tfa.pl
#    Erich Schwarz, 11/25/01
# 
# Purpose: Large-scale homolog searches can be aided by batch Entrez,
#    but neither FASTA nor GenBank downloads really are ideal.
#    The FASTA files give gi|N but not real gene names; they also
#      can have vast comment lines.  The GenPept files have all
#      data I want but not in FASTA format.  Thus, this script is
#      designed to take a GenPept file and convert it into a FASTA
#      file with the following characteristics: the ">" line has
#      the *gene name* as its leading string; the gi|N identifier
#      is given later; and a clean FASTA sequence follows.

# 1a. If not given a genbank file as argument, ask for its name.

if ($#ARGV != 0) {
	print "Required: input GenBank file!\n";
	print "What will input file be? ";
	$infile = <STDIN>;
} else {
	$infile = $ARGV[0];
}
chomp ($infile);
$outfile = ($infile . ".tfa");
print "The input file is $infile; the output file $outfile\n";

# 1b. This is meant to allow outputs
#     in either of two formats (sometimes one's better than 
#     the other, sometimes the reverse).

# Initialize a control variable which I'll use to prevent
# text from being scanned as coding sequences.

$nt_read_toggle = 0;

# 2. Open the genpept and future .tfa files or die.

open (INFILE, "$infile") || die "GenBank file $infile not found. $!\n";
open (OUTFILE, ">$outfile") || die "Couldn't open file $outname. $!\n";

# 3. Extract everything I want in a format which is of real use.

while (<INFILE>) {
	if ($_ =~ /LOCUS/) {
		$nt_read_toggle = 0;
		}

# get number out of lines like:
# VERSION     BF264750.2  GI:13261709

	elsif ($_ =~ /^VERSION.+GI\:([0-9]+)$/) 
        {
		$nt_read_toggle = 0;
		$gi_number = $1;
		chomp($gi_number);
		$gi_number =~ s/[^0-9]//g;

                $dev_stage_name = "unspecified";
                $tissue_type_name = "unspecified";
                $sex_name = "unspecified";
                $clone_lib_name = "unspecified";
                $species = "no_species_name";

                $reading_clone_lib = "no";
	}
	elsif ($_ =~ /ORGANISM[\s]+[.]*/) {
		$species = $_;
		chomp($species);
		$species =~ s/ORGANISM[\s]+//g;
	}

# extract data from lines like:
# /dev_stage="Parasitic Adults"
# or
# /tissue_type="pooled"
# or
# /clone_lib="unpublished oligo-capped cDNA library, stage
# L2"
# or
# /sex="Hermaphrodite"

        elsif ($_ =~ /dev_stage=/)
        {
                $dev_stage_name = $_;  
                chomp($dev_stage_name);
                $dev_stage_name =~ s/[\W]*dev_stage=//;
                $dev_stage_name =~ s/\"//g;
        }

        elsif ($_ =~ /tissue_type=/) 
        {   
                $tissue_type_name = $_;  
                chomp($tissue_type_name);
                $tissue_type_name =~ s/[\W]*tissue_type=//;
                $tissue_type_name =~ s/\"//g;
        }
        elsif ($_ =~ /sex=/)
        {
                $sex_name = $_;
                chomp($sex_name);
                $sex_name =~ s/[\W]*sex=//;
                $sex_name =~ s/\"//g;
        }

#        elsif ($_ =~ /clone_lib=/)
#        {
#                $clone_lib_name = $_;
#                chomp($clone_lib_name);
#                $clone_lib_name =~ s/[\W]*clone_lib=//;
#                $clone_lib_name =~ s/\"//g;
#                $reading_clone_lib = "yes";
#        }
#        elsif (($_ =~ /[^=]/) && ($reading_clone_lib eq "yes"))
#        {
#                $clone_lib_addon = $_;
#                chomp($clone_lib_addon);
#                $clone_lib_addon =~ s/[\W]*clone_lib=//;
#                $clone_lib_addon =~ s/\"//g;
#                $clone_lib_name = $clone_lib_name . " " . $clone_lib_addon;
#        }
#        elsif (($_ =~ /=/) || ($_ =~ /BASE COUNT/))
#        {
#                $reading_clone_lib = "no";
#        }

	elsif ($_ =~ /ORIGIN/) 
        {
            print OUTFILE ">";
            print OUTFILE $gi_number;
            print OUTFILE "   ";
            print OUTFILE $species;
            print OUTFILE "  ";

            unless ($dev_stage_name eq "unspecified" || $dev_stage_name eq "") 
            {
                print OUTFILE "    Dev. stage: ";
                print OUTFILE $dev_stage_name;
            }
            unless ($tissue_type_name eq "unspecified" || $tissue_type_name eq "")
            {
                print OUTFILE "    Tissue type: ";
                print OUTFILE $tissue_type_name;
            }
            unless ($clone_lib_name eq "unspecified" || $clone_lib_name eq "") 
            {
                print OUTFILE "    Clone lib.: ";
                print OUTFILE $clone_lib_name;
            }
            unless ($sex_name eq "unspecified" || $sex_name eq "") 
            {
                print OUTFILE "    Sex: ";
                print OUTFILE $sex_name;
            }
            print OUTFILE "\n";
            $nt_read_toggle = 1;
	}

	elsif ($_ =~ /[a-zA-Z]+/) 
        {
	    if ($nt_read_toggle == 1) 
            {
	        $dna_seq_line = $_;
	        chomp($dna_seq_line);
	        $dna_seq_line =~ s/[\d]*//g;
                $dna_seq_line =~ s/[\s]*//g;
                $dna_seq_line =~ s/[\W]*//g;
		$dna_seq_line =~ tr/a-z/A-Z/;
		print OUTFILE ($dna_seq_line . "\n");
	    }
	}

        elsif ($_ =~ /\\/) 
        {
            if ($nt_read_toggle == 1) 
            {  
                $nt_read_toggle = 0;
            }
        }

}
