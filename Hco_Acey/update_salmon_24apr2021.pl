#!/usr/bin/env perl

use strict;
use warnings;
use autodie;

while (my $input = <>) {
    chomp $input;

    $input =~ s#cd /pylon5/mc5phvp/schwarze/Acey/2019.05/salmon#cd \$PROJECT/Acey/2021.04.22/salmon#;
    $input =~ s#\. \$SCRATCH/anaconda2/etc/profile.d/conda.sh#\. \$PROJECT/anaconda3/etc/profile.d/conda.sh#;
    $input =~ s#\$SCRATCH/Acey/2019.05/salmon/dbs/ancylostoma_ceylanicum.PRJNA231479.WBPS13.CDS_transcripts#\$PROJECT/Acey/2021.04.22/salmon/dbs/Acey_WBPS14_gentrome_index#;
    $input =~ s#\$SCRATCH/Acey/2019.05/salmon/annots/Acey_PRJNA231479.WBPS13.cds2gene.tsv.txt#\$PROJECT/Acey/2021.04.22/annots/Acey_PRJNA231479.WBPS14.cds2gene.tsv.txt#;

    $input =~ s#/pylon5/mc5phvp/schwarze#\$PROJECT#g;

    $input =~ s#salmon_0.13.1#salmon_1.3.0#;
    $input =~ s#--validateMappings#--seqBias --gcBias --validateMappings#g;
    $input =~ s#--seqBias --gcBias --geneMap#--geneMap#g;

    $input =~ s#2019.05.28#2021.04.24#g;
    $input =~ s#\$SCRATCH#\$PROJECT#g;
    $input =~ s# --numBootstraps \d+##;

    print "$input\n";
}
