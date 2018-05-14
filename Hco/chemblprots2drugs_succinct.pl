#!/usr/bin/env perl

# chemblprots2drugs_succinct.pl -- Erich Schwarz <emsch@its.caltech.edu>, 3/5/2011.
# Purpose: link data from four different MySQL outputs into a single *succinct* table linking ChEMBL protein names to intelligible drug targets.

use strict;
use warnings;
use Getopt::Long;

my $data_ref;

my $drug_chembl_id    = q{};   
my $drug_common_name  = q{}; 
my $drug_formula_name = q{};
my $drug_number       = q{};
my $prot_chembl_id    = q{};
my $prot_number       = q{};
my $confidence_score  = q{};
my $conf_threshold    = 0;

my %opts = ();

GetOptions(
    "wanted=s" => \$opts{'wanted'},
    "t1=s"     => \$opts{'table_1'},
    "t2=s"     => \$opts{'table_2'},
    "t3=s"     => \$opts{'table_3'},
    "t4=s"     => \$opts{'table_4'},
    "c=i"      => \$conf_threshold,
    "help"     => \$opts{'help'},
);

if (    $opts{'help'} 
     or (! $opts{'wanted'}  )
     or (! $opts{'table_1'} ) 
     or (! $opts{'table_2'} )
     or (! $opts{'table_3'} )
     or (! $opts{'table_4'} ) 
     or ( $conf_threshold != int($conf_threshold) ) 
     or ( $conf_threshold < 0 ) 
     or ( $conf_threshold > 9 )
   ) { 
    die "Format: chemblprots2drugs_succinct.pl",
        " --wanted [file with lines starting with wanted proteins]",
        " --t1 [Table 1] --t2 --t3 --t4 [other tables 2-4] -h|--help",
        " -c [optional confidence threshold, 1-9]",
        "\n",
        ;
}

if (! $conf_threshold) { 
    $conf_threshold = 1;
}

open my $WANTED, '<', $opts{'wanted'} or die "Can't open file with list of wanted ChEMBL prots.\n";
while (my $input = <$WANTED>) {
    chomp $input;
    if ( $input =~ /\A (CHEMBL\d+) \s /xms ) { 
        $prot_chembl_id = $1;
        $data_ref->{'wanted_prot_CHEMBL_id'}->{$prot_chembl_id} = 1;
    }
}
close $WANTED or die "Can't close filehandle to list of wanted ChEMBL prots.\n";

open my $TAB4, '<', $opts{'table_4'} or die "Can't open table 4, $opts{'table_4'}.\n";
while (my $input = <$TAB4>) { 
    chomp $input;

# Sample input:
#    | tid    | chembl_id     |
#    | 102973 | CHEMBL1075023 |

    if ( $input =~ / \A 
                     (?: \| \s+ ) 
                     (\d+) 
                     (?: \s+ \| \s+ ) 
                     (CHEMBL\d+) 
                     (?: \s+ \| ) 
                   /xms ) { 
        $prot_number    = $1;
        $prot_chembl_id = $2;

        # Should be unique mappings of different ID types.

        # This is not general-purpose -- it's censoring most of the data --
        #     but, where we want only to know about a few proteins, it does cut 
        #     way down on unwanted lines of final output.

        if ( exists $data_ref->{'wanted_prot_CHEMBL_id'}->{$prot_chembl_id} ) { 
            $data_ref->{'prot_number'}->{$prot_number}->{'prot_CHEMBL_id'}    = $prot_chembl_id;
            $data_ref->{'wanted_prot_number'}->{$prot_number} = 1
        }
    }
}
close $TAB4 or die "Can't close filehandle to table 4, $opts{'table_4'}.\n";

open my $TAB3, '<', $opts{'table_3'} or die "Can't open table 3, $opts{'table_3'}.\n";
while (my $input = <$TAB3>) {
    chomp $input;

# Sample input:
#    | molregno | confidence_score | tid    |
#    |   106865 |                0 |  22224 |

    if ( $input =~ / \A 
                     (?: \| \s+ ) 
                     (\d+) 
                     (?: \s+ \| \s+ ) 
                     (\d+) 
                     (?: \s+ \| \s+ ) 
                     (\d+) 
                     (?: \s+ \| ) 
                   /xms ) { 
        $drug_number      = $1;
        $confidence_score = $2;
        $prot_number      = $3;

        # Not necessarily unique mappings; there can be many drugs/protein, and vice versa.
        # Don't bother with confidence scores of '0'!  $conf_threshold default value is '1'; can go up to '9'.
        # And don't enter these data if 'prot_number' wasn't entered, just above.
        #    (Which, in turn, should only happen if $prot_chembl_id was in our wanted file.)

        if ( ( $confidence_score >= $conf_threshold ) and ( exists $data_ref->{'wanted_prot_number'}->{$prot_number} ) ) { 
            $data_ref->{'prot_number'}->{$prot_number}->{'drug_number'}->{$drug_number} = 1;
            $data_ref->{'wanted_drug_number'}->{$drug_number} = 1;
        }
    }
}
close $TAB3 or die "Can't close filehandle to table 3, $opts{'table_3'}.\n";

open my $TAB1, '<', $opts{'table_1'} or die "Can't open table 1, $opts{'table_1'}.\n";
while (my $input = <$TAB1>) {
    chomp $input;
                     
# Sample input:
#     | molregno | chembl_id     |
#     |     2223 | CHEMBL41      |
        
    if ( $input =~ / \A
                     (?: \| \s+ )
                     (\d+)
                     (?: \s+ \| \s+ )
                     (CHEMBL\d+)
                     (?: \s+ \| )
                   /xms ) {
        $drug_number    = $1;
        $drug_chembl_id = $2;
        
        # Should be a unique mapping, both ways.
        if ( exists $data_ref->{'wanted_drug_number'}->{$drug_number} ) {
            $data_ref->{'drug_number'}->{$drug_number}->{'drug_chembl_id'}    = $drug_chembl_id;
            $data_ref->{'wanted_drug_CHEMBL_id'}->{$drug_chembl_id} = 1;
        }
    }
}
close $TAB1 or die "Can't close filehandle to table 1, $opts{'table_1'}.\n";

open my $TAB2, '<', $opts{'table_2'} or die "Can't open table 2, $opts{'table_2'}.\n";
while (my $input = <$TAB2>) {
    chomp $input;
    
# Sample input:
#    | chembl_id | pref_name | compound_name |
#    | CHEMBL6329    | NULL       | 2-[4-(2-Chloro-benzoyl)-3-methyl-phenyl]-2H-[1,2,4]triazine-3,5-dione |
#    | CHEMBL41      | FLUOXETINE | Methyl-[3-phenyl-3-(4-trifluoromethyl-phenoxy)-propyl]-amine          | 

    if ( $input =~ / \A
                     (?: \| \s+ ) 
                     (CHEMBL\d+)
                     (?: \s+ \| \s+ )
                     ([A-Z]+)
                     (?: \s+ \| \s+ )
                     (\S.+\S)
                     (?: \s+ \| )
                   /xms ) { 
        $drug_chembl_id    = $1;
        $drug_common_name  = $2;
        $drug_formula_name = $3;

        if ( exists $data_ref->{'wanted_drug_CHEMBL_id'}->{$drug_chembl_id} ) { 
            # Avoid GREAT RUNES in common names.
            $drug_common_name  =~ tr/[A-Z]/[a-z]/;

            # Should be a unique mapping.
            $data_ref->{'drug_chembl_id'}->{$drug_chembl_id}->{'drug_common_name'} = $drug_common_name;

            # *Not* a unique mapping.
            # Also, do not enter chemical names redundant with the common name; make cross-check case-insensitive ('/i').
            if ( $drug_formula_name !~ /\A$drug_common_name\z/i ) { 
                $data_ref->{'drug_chembl_id'}->{$drug_chembl_id}->{'drug_formula_name'}->{$drug_formula_name} = 1;
            }
        }
    }
}
close $TAB2 or die "Can't close filehandle to table 2, $opts{'table_2'}.\n";

my @prot_nums1 = grep { $_ =~ /\A\d+\z/ } sort keys %{ $data_ref->{'prot_number'} };
foreach my $prot_num1 ( @prot_nums1 ) { 
    $prot_chembl_id = $data_ref->{'prot_number'}->{$prot_num1}->{'prot_CHEMBL_id'};
    my @drug_numbers = sort keys %{ $data_ref->{'prot_number'}->{$prot_num1}->{'drug_number'} };
    foreach my $drug_num1 (@drug_numbers) {
        # Gather data: 
        if (! exists $data_ref->{'drug_number'}->{$drug_num1}->{'drug_chembl_id'} ) { 
            die "Inexplicable failure to get drug ChEMBL ID number.\n";
        }
        $drug_chembl_id   = $data_ref->{'drug_number'}->{$drug_num1}->{'drug_chembl_id'};
        $drug_common_name = q{};
        if ( exists $data_ref->{'drug_chembl_id'}->{$drug_chembl_id}->{'drug_common_name'} ) {
            $drug_common_name = $data_ref->{'drug_chembl_id'}->{$drug_chembl_id}->{'drug_common_name'};
        }
        my @drug_formula_names    = sort keys %{ $data_ref->{'drug_chembl_id'}->{$drug_chembl_id}->{'drug_formula_name'} };
        if ( $drug_common_name ne 'null' ) {
            push @drug_formula_names, $drug_common_name;
            @drug_formula_names = sort @drug_formula_names;
        }
        my $drug_formula_namelist = join '; ', @drug_formula_names;
        $drug_formula_namelist    = q{"} . $drug_formula_namelist . q{"};
        $data_ref->{'final_prot_name'}->{$prot_chembl_id}->{'final_drug'}->{$drug_formula_namelist} = 1;
    }
}
                                  
foreach my $final_prot ( sort keys %{ $data_ref->{'final_prot_name'} } ) { 
    my @final_drugs = sort keys %{ $data_ref->{'final_prot_name'}->{$final_prot}->{'final_drug'} };
    my $final_drug_text = join '; ', @final_drugs;
    print "$final_prot\t$final_drug_text\n";
}

