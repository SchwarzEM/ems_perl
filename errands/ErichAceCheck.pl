#!/usr/bin/perl -w

# ErichAceCheck.pl
# by Wen Chen, 2/2002
# some minor tweaks by Erich Schwarz

#-------------------Type out the purpose of the script-------------
print "This program checks Locus and Sequence objects.\n";
print "This is designed for Erich's .ace files for gene annotation.\n";
print "-*-------------*--------------*--------------*----------\n\n";


#------------------Make dictionaries------------------------------
print "Make Locus dictionary ...";
@locus = MakeDictionary ("/home/schwarz/dict/l_Locus.ace");
print "Done.\n";
$lo_count=@locus;

print "Make Sequence dictionary ...";
@sequence = MakeDictionary ("/home/schwarz/dict/l_Sequence.ace");
print "Done.\n";
$se_count=@sequence;

print "Make Strange Sequence dictionary ...";
@strange_seq = MakeDictionary ("/home/schwarz/dict/l_StrangeSeq.ace");
print "Done.\n";
$stse_count=@strange_seq;

print "Make Gene_class dictionary ...";
@gene_class = MakeDictionary ("/home/schwarz/dict/l_Gene_class.ace");
print "Done.\n";
$ge_count=@gene_class;

#-----------Dictionaries made, now check .ace files--------------
print "What .ace file do you want to check? ";
chomp($Ace_name=<stdin>);
print "What do you want to call the output file? ";
chomp($Out_name=<stdin>);

print "Start checking ... \n\n";
open (IN, "$Ace_name") || die "can't open $!";	
open (OUT, ">$Out_name") || die "can't open $!";	
$Line=<IN>;
$l=1;
while ($Line=<IN>) {
        $l++;
	chomp ($Line);
	if ($Line ne "") {
           @word = split ('"', $Line);
           $length = @word;
           if ($length >= 2) {
	      @w = split (' ', $word[0]);
	      $word[0] = $w[0];
	      $word[1] = "\"$word[1]\"";
	      $step=0; 
	      $match=0;
	      if (($word[0] eq "Locus") || 
		  ($word[0] eq "Other_name") || 
		  ($word[0] eq "Old_name"))  {	
		 while (($step < $lo_count) && ($match == 0)) {
		    if ($word[1] eq $locus[$step]) {
			$step=$lo_count+1;
			$match = 1; #found the locus in dictionary
		    } else {
			$step++;
		    }
		 }
	      } elsif (($word[0] eq "Sequence") || 
		       ($word[0] eq "Genomic_sequence")) {
		 while (($step < $se_count) && ($match == 0)){
		    if ($word[1] eq $sequence[$step]) {
			$step=$se_count+1;
			$match = 1; #found the sequence in valid dictionary
		    } else {
			$step++;
		    }
		 }
		 $step = 0;
		 while (($step < $stse_count) && ($match == 0)){
		    if ($word[1] eq $strange_seq[$step]) {
			$step=$stse_count+1;
			$match = 2; #found the sequence in invaid dictionary
		    } else {
			$step++;
		    }
		 }
	      } elsif (($word[0] eq "Gene_Class") || 
		       ($word[0] eq "Gene_class")) {
		 while (($step < $ge_count) && ($match == 0)){
		    if ($word[1] eq $gene_class[$step]) {
			$step=$ge_count+1;
			$match = 1;
		    } else {
			$step++;
		    }
		 }
       	      } else { $match = 1; }
	      if ($match == 0) {
		  print OUT "Line $l is new:  $Line\n"; 		
	      } elsif ($match == 2) {
		  print OUT "Line $l sequence might be invalid:  $Line\n";    
	      }	  
	  }  #End of if ($length >= 2) 
       }  # End of if {$Line ne "") 
}	
close (IN);
close (OUT);
print "\n";
print "Totally $l lines read. Thank you for using aceChecker!\n";


#------------------------sub routines-----------------------------
sub MakeDictionary {
	my $filename = shift(@_);
	my @Type = "null";
	print "file name : $filename \n";
	open (IN, $filename) || die "can't open $!";	
	$Line=<IN>;
	$l=1;
	while ($Line=<IN>) {
		chomp ($Line);
		@fields = split / : /, $Line; 
		$Type[$l-1] = $fields[1];
		$l++;
	}	
	close (IN);
	return @Type;
}
