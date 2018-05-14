#!/usr/bin/perl -w

# match_sequences_to_interpro.pl
# Erich Schwarz, 1/29/02

# Goal: map sequences, e.g., prion gene sequences, to 
#       Interpro text via Tablemaker file.

print "Input name of sequence list: ";
$sequence_list = <STDIN>;
chomp ($sequence_list);

print "Input name of sequence->Interpro Tablemaker file: ";
$seq2ipro_table = <STDIN>;
chomp ($seq2ipro_table);

print "Input Interpro term list: ";
$interproterms = <STDIN>;
chomp ($interproterms);

$output_file = $sequence_list . "seq2interpro_output";

system "touch $output_file";

open (OUTPUT_FILE, ">>$output_file")
    || die "Can't open output file $output_file. $! \n";

open (SEQLIST, "$sequence_list") 
    || die "Can't open sequence list $sequence_list. $! \n";

# Format of sequence list is (e.g.):

# W03D2.10
# W06B11.2
# Y56A3A.4

while (<SEQLIST>) 
{
    $sequence = $_;
    chomp ($sequence);
    foreach ($sequence) 
    {
        open (SEQ2IPRO, "$seq2ipro_table")
            || die "Can't open sequence->Interpro Tablemaker file $seq2ipro_table. $! \n";

# Format of sequence->Interpro Tablemaker file is (e.g.):

# "W03D2.10"	"WP:CE18320"	"INTERPRO:IPR003002"	
# "W06B11.2"	"WP:CE05044"	"INTERPRO:IPR001313"	"GO:0003723"
# "Y56A3A.4"	"WP:CE22575"	"INTERPRO:IPR000166"	"GO:0003677"

        while (<SEQ2IPRO>) 
        {
            $seq2ipro_line = $_;
            chomp ($seq2ipro_line);
            foreach ($seq2ipro_line) 
            {
                if ($seq2ipro_line =~ /^\"$sequence\".+\"INTERPRO.(IPR\n+)\"/) 
                {
                    $ipro_gene_tag = $1;
                    open (IPROTERMS, "$interproterms")
                        || die "Can't open Interpro term list $interproterms. $! \n";

# Format of Interpro term list is (e.g.):

# IPR003002 7TM chemoreceptor subfamily 1
# IPR000166 Histone-fold/TFIID-TAF/NF-Y domain
# IPR001313 Pumilio-family RNA binding domains (aka PUM-HD, Pumilio homology domain)

                    while (<IPROTERMS>) 
                    {
                        $iproterm_line = $_;
                        chomp ($iproterm_line);
                        foreach ($iproterm_line) 
                        {
                            $correctly_found_match = "no_match_found";
                            if ($iproterm_line =~/^$ipro_gene_tag/) 
                            { 
                                print OUTPUT_FILE "$sequence\t$iproterm_line\n";
                                $correctly_found_match = "match_was_found";
                            }
                        }
                        if ($correctly_found_match eq "no_match_found") 
                        {
                            die "Failed at Interpro term list acquisition step. $! \n";
                        }
                    close IPROTERMS;
                    }
                }
            }
        }
        close SEQ2IPRO;
    }
}

close SEQLIST;
close OUTPUT_FILE;
