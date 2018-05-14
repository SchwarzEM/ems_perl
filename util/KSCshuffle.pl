#! /usr/bin/perl
#
# KSCshuffle.pl
# KSCschuffle contains a Pseudo Random Number Generator (KSC_PRNG) that 
# generates pseudo random Bytes and fractions (6 Byte resolution) and
# shuffles lists.
# 
# The KSC_PRNG takes a Byte string as seed and can exploit any "randomness" 
# in the seed string. That is, to get a seed equivalent to a 32 Byte 
# (256 bit) sequence with a uniform probability density, you can just
# enter an ASCII sequence of 45 cards symbols drawn from a full deck
# of playing cards (with reshuffling between draws). Or you can use the
# recorded noise of you computer fan.
#
# The algorithm used is completely idosyncratic and the result of a hobby 
# of mine. It is untested and I have found no descriptions of it anywhere. 
# My mathematical skills are very limited so don't count on this to be safe
# in ANY way. If you still haven't understood it, THIS IS AN UNPROVEN PRNG,
# DON'T EXPECT IT TO BE "GOOD" IN ANY SENSE.
#
# That said, my understanding until now is that it has a large cycle time
# (O(2**4100)) and the output seems to have a rather uniform distribution.
# Anyway, it will undoubtedly be better (and slower) than the braindead 
# PRNG in the standard distribution of perl.
# 
# How to use it.
# Command Line:
# cat File | KSCshuffle.pl <seed> > resultFile
# 
# Inside Perl scripts:
# require 'KSCshuffle.pl';
# $PRNG = new KSC_PRNG(<seed>);		# Seed PRNG
# $RandomByte = $PRNG->nextByte;	# Get a Pseudo Random Byte
# $RandomFract = $PRNG->fraction;	# Get a Pseudo Random fraction ([0, 1>)
# $PRNG->Shuffle(\@List);			# Shuffle a list
#
# If no seed is entered, the function RealNoise is called that attemps to
# extract some noisy data from your system. You are adviced to connect it
# to a microphone, webcam or some other device with unpredictable output.
#
# 08 Jan 2002 - IMPORTANT: Changed PRNG output to every THIRD byte. Using
#               every SECOND byte allows for attacks!
#
############################################################################
#                                                                          #
#   Copyright (c) 1999                                                     #
#	Rob van Son                                                        #
#	Rob.van.Son@hum.uva.nl                                             #
#	Rob.van.Son@workmail.com                                           #
#	Institute of Phonetic Sciences/IFOTT                               #
#	University of Amsterdam                                            #
#	Herengracht 338                                                    #
#	1016CG Amsterdam                                                   #
#	The Netherlands                                                    #
#                                                                          #
#   This program is free software; you can redistribute it and/or modify   #
#   it under the terms of the GNU General Public License as published by   #
#   the Free Software Foundation; either version 2 of the License, or      #
#   (at your option) any later version.                                    #
#                                                                          #
#   This program is distributed in the hope that it will be useful,        #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of         #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          #
#   GNU General Public License for more details.                           #
#                                                                          #
#   You should have received a copy of the GNU General Public License      #
#   along with this program; if not, write to the Free Software            #
#   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.              #
#   Even if we assume                                                      #
#                                                                          #
############################################################################
#
# KSC Pseudo-Random Number Generator
# 
# How does it work. 
# The basic operation is the scrambler (pseudo-code):
# NewValue[i] = MappingTable[OldValue[i] XOR NewValue[i-1]]
# with NewValue[-1] = NewValue[$#Values] (i.e., the last value)
# 
# The MappingTable is a (pseudo-) random remapping of values.
# 
# This operation is reversible and maximally dispersive. Flip one bit
# in the original string and after two rounds each bit has a 50/50 chance 
# of being different from the scrambled version of the original string.
# 
# The scrambler is used in two ways. First it is used to scramble the 
# plaintext so every (scrambled) bit depends on every bit of the 
# plaintext (this is a bit like the Fourrier and Walsh transforms). 
# Second, it is used to construct a pseudo-random number generator (PRNG). 
# The generator consists of a seed string that is repeatedly scrambled. 
# Each time a new mapping table is constructed from the scrambled string. 
# After 4 rounds, every seed string will lead to a flat bit distribution
# (except for a set of "regular" inputs).
# 
# The scrambler is completely reversible. This means that all information
# in the original string is retained. If the seed string of the PRNG is N 
# bytes long, the number of possible input strings is 2**(8*N). The 
# expected cycle time of the PRNG is of the order of 2**(4*N). Actually, 
# the cycle time will be determined by the length of the INTERNAL string 
# which will be at least 1024 bytes.
# 
# KSC takes string and seeds the PRNG with it. A 50 character random
# string with only lower-case characters and numbers (i.e., [a-z0-9])
# is equivalent to a uniformly distributed 256 bit seed. Note that the 
# 50 characters should be realy RANDOM. That is, Scrabble pieces drawn 
# repeatedly or a 36 hole roulette. Do NOT use plain English text. Even 
# a Jane Austen novel contains only 1 BIT of information for every 5 
# BYTES (and that is an optimistic estimate, realy, I calculated it). 
# 
# How to get at a "good" random 256 bit seed:
# 
# Toss a coin 256 times [HT]
# Throw a dice 99 times [1-6]
# Get 77 random digits  [0-9]
# Draw 70 cards from a full deck [1-9PQKA] (value only, reshuffle!)
# Draw 64 cards from a small deck [7-9PQKA][HDSC] (value&color, reshuffle!)
# Record the order of a shuffled deck of 58 DISTINCT cards (value&color)
# Draw 55 Scrabble "coins" from a bowl [a-z] (replace characters and mix)
# Get 50 random lowercase keys [0-9a-z]
# Draw 45 cards from a full deck [1-9PQKA][HDSC] (value&color, reshuffle!)
# Get 42 random lower- and uppercase keys [0-9a-z!@#$%^&*()A-Z]
# 
# Internaly the PRNG will use at least 1024 bytes. A 50 character PRNG seed 
# has 2**256 possible sequences. The number of distinct "states" of the 
# 1024 byte string is ~2**8192. The expected, average, cycle time of the 
# PRNG will be ~2**4096 with a very sharply peaked distribution. Only a 
# vanishing proportion of sequences will be either shorter than ~2**4095 
# or longer than ~2**4097. 2**4096 is also the expected number of distinct 
# repeating sequences. This means that the probabillity that two of the 
# 256 bit keys are part of the same sequence is P ~ 1/2**3048. Which can 
# be said to be "unlikely". 
# Although these numbers look impressive, they are meaningless unless it 
# can be made plausible that "short", non-uniform (i.e., compressible) 
# seed sequences will indeed distribute more or less evenly over the space 
# of possible sequences. In other words, the likely seed sequences will 
# all be rather close together (in terms of the Hamming distance). The 
# PRNG will be much less "safe" if they stay close together after the 
# scramblings (more on this later).
# 
# To disconnect the output from the internals of the PRNG only every 
# SECOND value is outputed. After the last value of the string has been 
# given, the internal string is scrambled two rounds and output starts 
# again at the first position. 
#
#  
# Mathematical background (very naive speculations)
#   
# Scrambler
# The scrambler has two components: 1 the Mapping table, 2 the String
# Currently, the scrambler works on 8 bit values (i.e., bytes). But it 
# can easily be reconfigured to work with 4 bit or 16 values.
# 
# 1 The Mapping table is constructed by "shuffling" the numbers 0-MapLength.
# Using an 8 bit, 256 byte table (starting with $MapTable = "0..255"):
# 
# for($i = 255; $i > 0; --$i)
# {
# 	my $j = vec($MapTableSeed, $i, 8);
# 	my $Value = vec($MapTable, $i, 8);
# 	vec($MapTable, $i, 8) = vec($MapTable, $j, 8);
# 	vec($MapTable, $j, 8) = $Value;
# };
#
# 2 The scrambler takes the string and performs the following operation:
# 
# my $PreviousByte = $Stringlength - 1;  # Byte -1 == last Byte
# for($j=0; $j < 2; ++$j) # > 1 loop
# {
# 	for($i=0; $i < $Stringlength; ++$i)
# 	{
# 		# New value is this byte and the previous one XORed (-1==last)
# 		my $Value = vec($String, $i, 8) ^ vec($String, $PreviousByte, 8);
# 		# The new value is mapped
# 		vec($String, $i, 8) = vec($MapTable, $Value, 8);
# 		$PreviousByte = $i;
# 	};
# };
#
# This reversible operation scrambles the whole string after only 2 rounds.
# It preserves "information". The following conjectures assume that the
# mapping table is incompressible, i.e., cannot be replaced by a shorter
# function (e.g., Identity(), NOT(), MIRROR_BITS() or non-mixing bit positions).
# The scrambler has the following features:
# a If the string is incompressible, i.e., cannot be produced by a program
#   shorter than itself, the scrambled string too is incompressible in all 
#   but a vanishingly small proportion of cases. This follows from the 
#   reversibility of the operation. This means that a random string of bits 
#   is mapped to another random string of bits.
# b The average number of double rounds before the scrambler reproduces 
#   the same string again (i.e., the cycle time) is
#   2**($Stringlength*8/2), i.e., quite large. The average number of such 
#   cycles available too is 2**($Stringlength*8/2).
# c In spite of the above, the use of the eXclusive OR operation ensures
#   that there are some "special" sequences that do not "randomize". 
#   First, a string filled with all identical $MapTable[0] bytes (i.e. the 
#   byte that maps to 0) will be left unchanged. Second, there generally is 
#   a single pair of Bytes that maps onto itself. Third, there will be 
#   sequences that show "regular" behavior. These kinds of problem 
#   can be alleviated by using two different scramblers in series, instead 
#   of one.
# d If the string is long and itself periodic (period T bytes), the
#   scrambled string will (on average) have a period of length($MapTable)*T
#   bytes. (this is a very, very naive conjecture, based on 
#   sqrt(length($MapTable)) for each round). Note that this is the same 
#   kind of period doubling you see in the onset of deterministic chaos 
#   (c.f., fractals).
# e If we extend this naive reasoning somewhat further, we can conjecture
#   that if the original string was compressible to a size T, the scrambled
#   string will be compressible to a size length($MapTable)*T, unless the
#   string is descrambled first.
# f We can use point (d) to estimate that we need only  
#   log(length($Stringlength))/log(length($MapTable))
#   (double) scramble rounds to construct a string with a uniform 
#   distribution of bits and bytes. For a 256 byte table and a 1024 byte
#   string this corresponds to less than two calls to the scrambler
#   (actually, three loops of $j). However, ad c ensures that there is 
#   and almost trivial algorithm to map out the table if only two loops are
#   used.
#   
# Pseudo-Random Number Generator (PRNG)
#
# The PRNG used is just a scrambler, for which both the table and the 
# string are derived from a single Seed string. The PRNG is designed to 
# exploit any "randomness" available in the seed. A list of 45 two byte 
# ASCII playing card symbols (recorded from actually drawing cards) is 
# equally well used as a short 32 byte long random bitpattern drawn from a
# uniform distribution
# . 
# Because the Seed might be atypical (i.e., compressible), care must be
# taken that the mapping table is NOT compressible (e.g., non-randomizing).
# This is done by repeatedly calling the Randomize function. This function
# takes the string, uses it's first length($MapTable) bytes to construct
# a new Mapping table, and then scrambles the string. After four such
# "Randomizations" we can be reasonably sure that both the Mapping table
# and the String have uniform bit distributions (i.e., "incompressible"
# except by recreating the original seed).
# 
# A good PRNG is a program (algorithm) that expands a Seed string of bits
# to a longer bit-stream that is indistinguishable from a random bit-
# stream with uniform density. Indistinguishable is defined as the 
# impossibility to find an algorithm that has a reasonable chance of 
# distinguishing the PRNG sequence from a purely random sequence within 
# given "cost" bounds. That is, for all "figures of merit", d, there 
# should be NO algorithm that has a non-vanishing chance (i.e., 
# Prob > 1/P(d)) of finding a non-random aspect of the PRNG sequence 
# within D < P(d) steps. In this, P(d) is any polynomial function of d 
# (e.g., D ~ a0+a1*d+a2*d**2+...+ak*d**k) and d some increasing 
# function of the length of the Seed (e.g., it's size in bits). Such an 
# algorithm, if it would exist, is called efficient. In practise, we 
# would like that the probability that the PRNG sequence is "cracked" 
# at least decreases as ~1/2**SeedLength for a given cost or the cost 
# increases as fast or faster than ~2**SeedLength. In reality, such a 
# proof can hardly, if ever, be given. The most one can do is 
# specifying that a PRNG is good under some specific classes of attacks.
# 
# Our scrambler based PRNG generates distinct sequences for distinct Seeds.
# There are more than ~2**8192 possible internal states. Cycle times
# are expected to run in the ~2**4096 range for a total of ~2**4096 such
# cyclic sequences. As the state sequence is completely reversible back to 
# the original Seed string, we can be quite certain that each distinct seed 
# leads to a distinct output. Theoretically, the whole PRNG can be 
# reconstructed, and its future output predicted, with only 1068 output 
# bytes (NOTE: the PRNG outputs only every SECONG byte of the SECOND 
# scrambling). However, this would require comparing this output with ALL 
# 2**4096 possible sequences (or with all possible seeds, which is a
# smaller number), which  amounts to a specification of a GOOD PRNG.
# 
# Another point is how to PREDICT the output of the PRNG given its previous
# output. This point is actually connected to the previous one, that of
# distinguishing the PRNG sequence from random noise. One way to tackle
# this question is to link the scramble operation to a system that is known 
# to be "impossible" to compress or predict. 
# A Universal Turing machine (UTM) is such a system. As it is the archetype
# of a computer, it can be simulated by another one, but its calculation 
# cannot be "compressed" or reconstructed from fragments (at least, the 
# prospects of doing so are moot in general). If it could be proven that 
# the scrambling operation is equivalent to the operation of a UTM, then 
# the safety of the PRNG would be bolstered as less than half of its 
# internal activity is accessible for observation. In terms of UTM's, the 
# scrambler can be simulated by a finite state Turing machine with the 
# following restrictions:
# a A finite, circular tape.
# b Completely reversible operation.
# c Unidirectional, single step processing of the tape.
# d The machine state is identical to the previous tape symbol.
# 
# For instance, it is possible to construct mapping tables that allow the 
# Scrambler to execute reversible algorithms consisting of eXclusive Or, 
# Not, and bit-cycling operations, programmed in the string. 
# Together, all these features still do not constitute any proof of 
# "UTM equivalence" of the scrambler. Given the limited length of the 
# string on which everything is "calculated", it even has more of 
# something described by a Context-Sensitive Grammar. But there are 
# no efficient algorithms for general CS parsers either (at least, that is
# how I understood it).
# All these are only hints that the PRNG as defined here might be "good".
#
# 
# Attempt to some more rigorous foundations
# 
# As our standard of quality we choose a seed that is as unpredictable as
# a perfectly uniformly distributed 256 bit seed.
# 
# Guessing the future output of the scrambler should be at least as unlikely 
# as guessing a 256 bit uniformly distributed seed. 
# A perfectly good mixer of bits can be defined as:
# Sourcetext (s) <--MX--> Mixtext (t)
# Such that for every pair of functions A=a(s) and B=b(t) and scrambler X 
# it holds that:
# P(B|AX) == P(B|X) and P(A|BX) == P(A|X)
# That is, no information on 't' can be extracted from partial information
# on 's', and vice versa (independence). 
#	 
# There is no reason to assume that such perfect reversible scramblers 
# indeed exist. However, we can settle for less. We consider a scrambler 
# a good mixer with some figure of merrit 'm' if:
# P(B|AX) <= MAX(P(B|X), m) and P(A|BX) <= MAX(P(A|X), m)
# for most values (A,B) and functions (a,b). 'Most' meaning that
# the probability of hitting on values (A,B) and functions (a,b)
# that violate this rule WITHOUT knowing X is smaller than 'm'.
#    
# I don't know whether such good scramblers exist either, but I 
# will use a qualified version. A scrambler is a good mixer in 
# a' if the above holds at least for A'=a(s), B'=a'(t), and 
# m = P(sX) = P(Seed) >= 2**-256 (note: a'==b'). That is, knowing A' 
# gives you no useful knowledge of B'. 
# Note that in this case, it is not possible to obtain (partial) 
# information about an unknown scrambler itself when only A' and B' are 
# known.
# 
# It is required that the probability of obtaining a scrambler that 
# is NOT a good mixer in the choosen function a' is less than m. 
# 
# A non-good mixer can arise when not all bit-positions within Bytes
# mix. Our scrambler mixes in two directions. Horizontaly between Bytes 
# and vertically the bits within Bytes. The horizontal mixing uses the 
# eXclusive OR operation which is perfectly mixing. The vertical mixing 
# uses a random table lookup. Horizontally, our scrambler is a good 
# mixer if the number of mixing rounds is as large as the number of Bytes 
# in the input string. However, it is a good, horizontal, mixer of the 
# ODD Bytes after only a few rounds (around six rounds I think). Therefore 
# we use the interspersed PRN Bytes (i.e., a'(s) = ODD-Bytes-of(s)).
# 
# There are 256! possible scramblers (~2**1679). A seed selects one from an 
# (random) subset of 2**256 scramblers. The probability that a scrambler 
# does not mix all bits in a Byte is (for 1 non-mixing bit)
# P(1 bit non-mixing) = 8 . (128!)**2/256! ~ 2**-250. 
# "Higher order" non-mixing has much lower probabilities. This is clearly 
# NOT 256 bit seed equivalent. The probability that this can exploited
# depends on the function a' that describes the visible information. Our
# selection of a'(s) = ODD-Bytes-of(s) for a string length of 1024 Bytes 
# ensures that this loop-hole cannot be exploited.
# 
# In summary, this PRNG seems promising.
#
############################################################################
#
# How to get good seeds? This is a classical problem. If pure, uniform
# distributions are needed, only expensive solutions are available. 
# However, our PRNG is able to exploit ANY randomness. If a particular
# 2 kB string of bytes is as unpredictable as a 256 bit uniform distributed
# value, these strings are completely equivalent for our PRNG. You can
# record the sound of your computer fan, an /s/ sound, capture pictures
# with an on-line digital camera etc.. All these digital sources are
# largely unpredictable and can be used as seeds.
# 
# #######################################################################
# 
# Some auxilliary functions. This is an independent package. It can be
# used without the rest of the file present.
# 
package KSCauxilliary;

# Factorial
sub Fact   # ($i) -> factorial of i
{
	my $i = shift;
	my $result = 1;

	while($i > 1)
	{
		$result *= $i;
		--$i;
	};
	return $result;
}

#
# Real, real noise for seeds, i.e., unpredictable bytes
sub RealNoise # ()-> String
{
	my $Noise;
	#
	#
	# For Unix machines
	# $Noise = `ps -elf`;
	# Or even better
	$Noise = `ps -elf | gzip`;
	$Noise = `ps -elx | gzip` unless $Noise;
	$Noise = "".$$.$^T.localtime() unless $Noise; # UNSAFE, for testing only
	#
	# SGI Indy commands to get Noisy input
	# Record some sounds (e.g., an /f/, /s/, or your computer fan)
	# Make sure that the microphone actually records something.
	# (this one is mono, 16 bit, 8kHz, ulaw, 0.5 s)
	# 	`recordaifc -n 1 -s 16 -r 8000 -c ulaw -t 0.5 ran-00000.aifc`;
	# 	$Noise = `cat ran-00000.aifc | gzip`; # Get the compressed sound
	# 	`rm ran-00000.aifc`; # REMOVE THE FILE
	# Shoot a jpeg picture
	# Make sure the camera does capture something.
	#`vidtomem -f ran; dmconvert -f jfif ran-00000.rgb ran-00000.jfif`;
	#$Noise = `cat ran-00000.jfif`; # Read the picture
	#`rm ran-00000.rgb ran-00000.jfif`; # REMOVE THE PICTURES
	#
	return $Noise;
}
#
# Decode URL encoded arguments
sub URLdecode   # (URL encoded input) -> string
{
  my $input = shift || return "";
  my $output = "";
  my $char;
  my $Value;
  # Convert all "+" to " "
  $input =~ s/\+/ /g;
  # Convert all hexadecimal codes (%FF) to their byte values
  while($input =~ /\%([0-9A-F]{2})/)
  {
	  $output .= $`.chr(hex($1));
	  $input = $';
  };
  $output .= $input;  # The remaining part of $input
  #
  $output;
};

# Encode arguments as URL codes.
sub URLencode   # (input) -> URL encoded string
{
  my $input = shift || return "";
  my $output = "";
  my $char;
  my $inputLength = length($input);
  my $Value;
  my $i = 0;
  for($i=0; $i<$inputLength; ++$i)
  {
     $char = vec($input, $i, 8);
          if($char =~ /\s/)
          {  $output .= "+";}
          elsif($char =~ /\w/)
          {  $output .= $char;}
          else
          {  $output .= sprintf("%%%2.2x", ord($char));
          };
  };
  $output;
};

#
# Base64 MIME encoding of binary data
# Static tables for encoding and decoding
my $Base64Table = 
	"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
# The inverse hash table for decoding
my %Base64InvTable =
('A',0,'B',1,'C',2,'D',3,'E',4,'F',5,'G',6,'H',7,'I',8,'J',9,'K',10,'L',11,
 'M',12,'N',13,'O',14,'P',15,'Q',16,'R',17,'S',18,'T',19,'U',20,'V',21,'W',22,'X',23,
 'Y',24,'Z',25,'a',26,'b',27,'c',28,'d',29,'e',30,'f',31,'g',32,'h',33,'i',34,'j',35,
 'k',36,'l',37,'m',38,'n',39,'o',40,'p',41,'q',42,'r',43,'s',44,'t',45,'u',46,'v',47,
 'w',48,'x',49,'y',50,'z',51,'0',52,'1',53,'2',54,'3',55,'4',56,'5',57,'6',58,'7',59,
 '8',60,'9',61,'+',62,'/',63);
 
# Encode argument string in Base64. Add padding unless the second argument
# is non-zero (NoPad != 0). Does not split lines with CRLF
sub Base64encodeString	# (input string [, NoPad]) -> Base64 encoded String (NO CRLF)
{
	my $InputString = shift || undef;
	my $NoPad = shift || undef;
	#
	return "" unless($InputString);
	#
	my $OutputLength = 4 * int(length($InputString) / 3);
	$OutputLength += 4 if length($InputString) % 3;
	my $PadLength = !$NoPad && length($InputString) % 3 ? 
						3-length($InputString) % 3 : 0;
	#
	my $OutputString = "\000" x $OutputLength;
	# Convert triplets to quartets
	my $i;
	for($i = 0; $i < $OutputLength/4; ++$i)
	{
		my $BaseQuartet = "\000\000\000\000";
		#
		my $bit;
		my $Offset = 0;
		for($bit=0; $bit < 12; ++$bit)
		{
			vec($BaseQuartet, $bit+$Offset/3, 2) = 
					vec($InputString, $i*12+$bit, 2);
			++$Offset;
		};
		#
		my $j;
		for($j=0; $j < 4; ++$j)
		{
			vec($OutputString, 4*$i+$j, 8) = 
				vec($Base64Table, vec($BaseQuartet, $j, 8), 8);
		};	
	};
	# Pad with '=' to complete quartet
	for($i = $OutputLength-$PadLength; $i < $OutputLength; ++$i)
	{
		vec($OutputString, $i, 8) = ord('=');
	};
	#
	return $OutputString;
}

# Decode argument string in Base64 to plain byte string. Does not take care of
# CRLF sequences but it does remove padding characters.
sub Base64decodeString	# (Base64 string) -> Original bytes
{
	my $InputString = shift || undef;
	#
	return "" unless($InputString);
	#
	my $InputLength = length($InputString);
	my $OutputLength = int(3/4 * $InputLength);
	#
	my $OutputString = "\000" x $OutputLength;
	# Decode quartets
	my $CodeChar = "\000";
	my $i;
	my $Pad = ($InputString =~ s/\=/\=/sg);
	for($i=0; $i < $InputLength; ++$i)
	{
		vec($InputString, $i, 8) = $Base64InvTable{chr(vec($InputString, $i, 8))};
	};
	
	# Convert quartets to triplets  
	for($i = 0; $i < $OutputLength/3; ++$i)
	{
		my $Triplet = "\000\000\000";
		my $bit;
		my $Offset = 0;
		for($bit=0; $bit < 12; ++$bit)
		{
			vec($Triplet, $bit, 2) = 
					vec($InputString, $i*16+$bit+$Offset/3, 2);
			++$Offset; 
		};
		for($j=0; $j < 3; ++$j)
		{
			vec($OutputString, 3*$i+$j, 8) = 
				vec($Triplet, $j, 8);
		};
	};
	# Remove padded bytes
	while($Pad--) { chop $OutputString;};
	#
	return $OutputString;
}

#
#
# Calculate the mean entropy of a String using BitLength sized units.
# Prints the mean and minimum entropy of the symbols, symbol pairs 
# (variable distance 1-32, BitLength <= 8) and symbol triples (1-32&1-8, 
# BitLength <= 4). 
# It also prints an approximation of the expected mean entropy of a sample
# of the provided size for a uniform distributed source. This approximation 
# uses a Taylor expansion of the LOG term in the Entropy formula. The even 
# order terms are evaluated as Chi-square variables (and moments).
#
sub Entropy  # ($String, $BitLength) prints entropies
{
	my $String = shift;
	my $BitLength = shift || 4;
	#
	#
	my ($i, $j, $k, $l) = (0, 0, 0, 0);
	die "Cannot handle $BitLength size units\n" unless $BitLength <= 16;
	my $Jmax = $BitLength <= 8 ? 16 : 0;
	my $Lmax = $BitLength <= 4 ? 8 : 0;
	#
	my @Bits = ();
	# Initialize frequency table
	for($j=0; $j<$Jmax; ++$j)
	{
		@{$Bits[$j]}=();
		#
		for($l=0; $l<$Lmax; ++$l)
		{
			@{$Bits[$j][$l]}=(0);
		};
	};
	#
	# Determine frequencies
	my $NoBytes = length($String)*8/$BitLength;

	my ($Byte0, $Byte1, $Byte2) = (0, 0, 0);
	#
	for($i=0; $i<$NoBytes; ++$i)
	{
		$Byte0 = vec($String, $i, $BitLength);
		++$Bits[0][0][$Byte0];
		#
		next unless $i+$Jmax < $NoBytes;
		for($j=1; $j<=$Jmax; ++$j)
		{
			$Byte1 = (vec($String, $i+$j, $BitLength) << $BitLength) + $Byte0;
			++$Bits[$j][0][$Byte1];
			#
			next unless $i+$j+$Lmax < $NoBytes;
			for($l=1; $l<=$Lmax; ++$l)
			{
				$Byte2 = (vec($String, $i+$j+$l, $BitLength) << 2*$BitLength) + $Byte1;
				++$Bits[$j][$l][$Byte2];
			};
		};
	};

	# Calculate Entropies
	for($j=0; $j<=$Jmax; ++$j)
	{
		for($l=0; $l<=$Lmax; ++$l)
		{
			my $k;
			my $OrderP = $j ? ($l ? 3 : 2) : 1; 
			my $Cats = 2**($BitLength*$OrderP); # The expected number of categories
			my $N = 0;
			my $H = 0;
			my $ListCats = scalar(@{$Bits[$j][$l]});
			for($k=0; $k < $ListCats; ++$k)
			{
				$N += $Bits[$j][$l][$k];
			};
			for($k=0; ($k < $ListCats) && $N; ++$k)
			{
				$H -= $Bits[$j][$l][$k] ? 
						($Bits[$j][$l][$k]/$N) * log($Bits[$j][$l][$k]/$N)
						:0;
			};
			$Bits[$j][$l][0] = $ListCats > 1 ? $H/log(2) : 0;
			$Bits[$j][$l][1] = $N;
			$Bits[$j][$l][2] = 0;
			# An estimation of the correction to the uniform density based on
			# the sample size (N and number of categories). Based on a Taylor
			# expansion of the log2(F+e), using only the even terms in e that map
			# to a Chi-square distribution (dof = $Cats-1, 
			# E(Chi-sqr**n) = 2**n * n! * dof). 
			# NOTE: THIS FORMULA WILL NOT CONVERGE!
			$N = 1 unless $N > 0;
			for($k=0; $k < 2; ++$k)
			{
				$Bits[$j][$l][2] += 
							$Cats**$k * 2**$k * Fact($k) * ($Cats-1.0)
							 / (($k+1)*($k+2) * $N**($k+1));
			};
			
			$Bits[$j][$l][2] /= log(2); 
		};
	};
	#
	my @Entropies = ();
	@{$Entropies[0]} = ();
	# 0 Order, plain entropies
	for($i=0; $i<3; ++$i)
	{
		$Entropies[0][$i] = $Bits[0][0][$i];
	};
	$Entropies[0][3] = "na";
	#
	# 1 order, pairs entropies
	my @Mean = (0, 0, 0);
	my $MinH = $BitLength*(2);
	for($j=1; $j<=$Jmax; ++$j)
	{
		for($i=0; $i<3; ++$i)
		{
			$Mean[$i]  += $Bits[$j][0][$i];
		};		
		$MinH    = $Bits[$j][0][0] if $Bits[$j][0][0] < $MinH;
	};
	for($i=0; $i<3 && $Jmax; ++$i)
	{
		$Entropies[1][$i] = $Mean[$i]/$Jmax;
	};		
	$Entropies[1][3] = $MinH;
	#
	# 2 order, triples entropies
	@Mean = (0, 0, 0);
	$MinH = $BitLength*(3);
	for($j=1; $j<=$Jmax; ++$j)
	{
		for($l=1; $l<=$Lmax; ++$l)
		{
			for($i=0; $i<3; ++$i)
			{
				$Mean[$i]  += $Bits[$j][$l][$i];
			};		
			$MinH    = $Bits[$j][$l][0] 
						if $Bits[$j][$l][0] < $MinH && $Bits[$j][$l][1];
		};
	};
	for($i=0; $i<3 && $Jmax; ++$i)
	{
		$Entropies[2][$i] = $Mean[$i]/$Jmax;
		$Entropies[2][$i] /= $Lmax if $Lmax;
	};		
	$Entropies[2][3] = $MinH;
	#
	for($j=0; $j < 3; ++$j)
	{
		print "$j order mean: $Entropies[$j][0] bits, N=$Entropies[$j][1], ",
		"(E=", $BitLength*($j+1)-$Entropies[$j][2], ") min= $Entropies[$j][3]\n";
	};
	print "\n";
}

# #######################################################################
# 
# Single byte scrambling
package Scrambler;

#
# Construct a scrambler
# 
sub new  # (\$String, \$MapTableSeedPtr) -> $self
{
	my $type = shift;
	my $self = {};
	$self->{'String'} = shift || undef;
	my $MapTableSeedPtr = shift || $self->{'String'};
	# Construct other parts
	$self->{'MapTable'} = ShuffleMapTable($MapTableSeedPtr);
	$self->{'InverseMapTable'} = InvertMapTable(\$self->{'MapTable'});
	$self->{'Length'} = length(${$self->{'String'}});
	#	
	return bless($self);
}

# A new string
sub EnterString  # (\$StringPtr) -> \$StringPtr
{
	my $self = shift;
	$self->{'String'} = shift || undef;
	$self->{'Length'} = length(${$self->{'String'}});
	return $self->{'String'};
};

# Reversibly scramble a string of bytes
sub Scramble # $self-> (\$StringPtr, \$MapTablePtr) -> \$ScrambledStringPtr
{
	my $self = shift;
	my $StringPtr = $self->{'String'};
	my $MapTablePtr = \$self->{'MapTable'};
	my $rounds = 3;   # Number of times to scramble
	#
    my $Stringlength = $self->{'Length'};
	return $StringPtr if $Stringlength <= 1; # Do NOT scramble single bytes
	#
    my ($i, $j) = (0, 0);
    my $PreviousByte = $Stringlength - 1;  # Byte -1 == last Byte
	#
    # Do a Cyclic XOR of each byte with its predecessor 
    # followed by a lookup in a Random SliceByte Map 
    for($j=0; $j < $rounds; ++$j) # > 1 loop
    {
      for($i = 0; $i <  $Stringlength; ++$i)
      {
         # New value is this byte and the previous one XORed (-1==last)
         my $Value = vec($$StringPtr, $i, 8) ^ 
                     vec($$StringPtr, $PreviousByte, 8);
         # The new value is mapped
         vec($$StringPtr, $i, 8) = vec($$MapTablePtr, $Value, 8);
         $PreviousByte = $i;
      };
   };
   #
   return $StringPtr;
}

# Reversibly descramble a string of bytes
sub DeScramble # $self-> (\$StringPtr, \$InverseMapTablePtr) -> \$ScrambledStringPtr
{
	my $self = shift;
	my $StringPtr = $self->{'String'};
	my $InverseMapTablePtr = \$self->{'InverseMapTable'};
	my $rounds = 3;   # Number of times to scramble
	#
    my $Stringlength = $self->{'Length'};
	return $StringPtr if $Stringlength <= 1; # Do NOT scramble single bytes
	#
	my ($i, $j);
	#
    # Do a lookup in an inverse Random Byte Map followed by a Cyclic XOR 
    # of each byte with its predecessor  
    for($j=0; $j < $rounds; ++$j) # > 1 loop
    {
      for($i = $Stringlength-1; $i >= 0; --$i)
      {
         my $PreviousByte = $i ? $i-1 : $Stringlength-1;
         # The new value is mapped
		 my $Value = vec($$StringPtr, $i, 8);
         $Value = vec($$InverseMapTablePtr, $Value, 8);
         # Old value is this byte and the previous one XORed (-1==last)
         vec($$StringPtr, $i, 8) = $Value ^ vec($$StringPtr, $PreviousByte, 8);
      };
   };
   #
   return $StringPtr;
}
#
# Construct a New Slice Byte Map from a string
sub ShuffleMapTable  # (\$MapTableSeedPtr) -> $MapTable;
{
	my $MapTableSeedPtr = shift || return undef;
	#
    my $MapTable = "\000" x 256; # 256 Bytes
	my $MapSize = 256;
    #
    my $i = 0;
    #
	# Construct a 1..256 byte list
	for($i=0; $i < $MapSize; ++$i)
	{
		vec($MapTable, $i, 8) = $i;
	};
	#
	# Shuffle the byte slice list
    for($i = ($MapSize-1); $i > 0; --$i)
    {
       my $j = vec($$MapTableSeedPtr, $i, 8);
       my $Value = vec($MapTable, $i, 8);
       vec($MapTable, $i, 8) = vec($MapTable, $j, 8);
	   vec($MapTable, $j, 8) = $Value;
    };
    #
    return $MapTable;
}
#
# Invert a Slice Byte Map and construct a new Inverted Byte Map
sub InvertMapTable # (\$MapTablePtr) -> $InvertedMapTable
{
	my $MapTablePtr = shift || return undef;
	#
	my $InvertedMapTable = "\000" x 256; # 256 Bytes
	my $MapSize = 256;
	#
    my $i = 0;
    #
	# Invert Map Table
	for($i=0; $i < $MapSize; ++$i)
	{
		my $j = vec($$MapTablePtr, $i, 8);
		vec($InvertedMapTable, $j, 8) = $i;
	}
	return $InvertedMapTable;
}
#
# Randomize a string by using it repeatedly as the seed of the Maptable
sub Randomize   # () -> \$String
{
	my $self = shift;
	my $rounds = 4;
	#
	#
	for($i = 0; $i < $rounds; ++$i)
	{
		my $StringPtr = $self->{'String'};
		# Store InverseMaps
		push(@{$self->{'Derandomize'}}, $self->{'MapTable'});
		$self->{'MapTable'} = ShuffleMapTable($StringPtr);
		$self->{'InverseMapTable'} = InvertMapTable(\$self->{'MapTable'});
		$self->Scramble;
	}
	#
	return $self->{'String'};
}

#
# DeRandomize a string by using the stored inverse maps
# this function is primarily used to ensure complete reversbility
# and testing
sub DeRandomize   # () -> \$String
{
	my $self = shift;
	my $rounds = 4;
	#
	#
	for($i = 0; $i < $rounds; ++$i)
	{
		$self->DeScramble;
		$self->{'MapTable'} = pop(@{$self->{'Derandomize'}});
		$self->{'InverseMapTable'} = InvertMapTable(\$self->{'MapTable'});
	}
	#
	return $self->{'String'};
}

# Test the integrity of the code
sub Test # () -> 1/0
{
	my $String = "a" x 65;
	my $Char = vec($String, 0, 8);
	my $Test = new Scrambler(\$String);
	$Test->Randomize;
	my $ScrambledString = $String;
	$Test->DeRandomize;
	my $i;
	my $n=0;
	for($i = 0; $i < 65; ++$i)
	{
		return 0 if vec($String, $i, 8) != $Char;
		++$n if vec($ScrambledString, $i, 8) == $Char;
	};
	return $n < 3 ? 1 : 0;
}

#
####################################################################
#
# The Pseudo-Random Number Generator
package KSC_PRNG;

@ISA = qw(Scrambler);
@KSC_PRNG::Inherit::ISA = @ISA;

sub new  # ($Seed) -> $self
{
	my $type = shift;
	my $Seed = shift;
	#
	# Without a seed, use process and timing information as a seed
	$Seed = KSCauxilliary::RealNoise unless($Seed);
	# Expand Seed if it is short
	$Seed .= ("\000" x (1024-length($Seed))) if length($Seed) < 1024;
    my $self = new Scrambler(\$Seed);

	$self->Randomize;
	# Construct other parts
	$self->{'Rounds'} = 0;
	$self->{'Pointer'} = 0;
	$self->{'CurrentIndex'} = 0;  # Sensitive to wordlength
	#	
	return bless($self);
}

#
# Get the current index in the expansion
sub CurrentIndex
{
	my $self = shift;
	#
	return $self->{'CurrentIndex'};
}

#
# Get the next index in the expansion
sub NextIndex
{
	my $self = shift;
	#
	++$self->{'CurrentIndex'};
	#
	return $self->{'CurrentIndex'};
}

#
# Get the previous index in the expansion
sub PreviousIndex
{
	my $self = shift;
	#
	--$self->{'CurrentIndex'};
	#
	return $self->{'CurrentIndex'};
}

# 
# Return the i-th pseudo random byte
sub indexByte  #($i) -> PRN i
{
	my $self = shift;
	my $index = shift || 0;
	$self->{'CurrentIndex'} = $index;
	# Only every THIRD byte is used
	$index *= 3;  
	my $ByteLength = $self->{'Length'};
	#
	my $rounds = int($index / $ByteLength);
	my $pointer = $index - ($rounds * $ByteLength);
	#
	# Scramble/DeScramble to get the correct version
	while($rounds > $self->{'Rounds'})
	{ 
		$self->Scramble;
		++$self->{'Rounds'};
    };
	while($rounds < $self->{'Rounds'})
	{ 
		$self->DeScramble;
		--$self->{'Rounds'};
	};
	# Update the position in the String
	$self->{'Pointer'} = $pointer;
	#
	# Result
	return vec(${$self->{'String'}}, $pointer, 8); 
}

# Random numbers between [0, 1> to 48 bit!
sub fraction  # () -> number between 0 and 1, uses next position
{
	my $self = shift;
	#
	my $ByteLength = 6;
	#
	my $i;
	my $result = 0;
	my $power = 256;
	for($i = 0; $i < $ByteLength; ++$i)
	{
	   $result += ($self->nextByte())/$power;
	   $power *= 256;
	};
	return $result;
}
#
# Next number, This is inline code because it is used VERY
# often.
sub nextByte
{
	my $self = shift;
	my $Pointer = $self->{'Pointer'} + 3;
	#
	++ $self->{'CurrentIndex'};
	#
	my $ByteLength = $self->{'Length'};
	# Scramble to get the correct version
	if($Pointer >= $ByteLength)
	{ 
		$self->Scramble;
		++$self->{'Rounds'};
		$Pointer = 0;	
    };
	$self->{'Pointer'} = $Pointer;
	#
	# Result
	return vec(${$self->{'String'}}, $Pointer, 8); 
}
#
# Previous number
sub previousByte
{
	my $self = shift;
	my $Pointer = $self->{'Pointer'} - 3;
	#
	-- $self->{'CurrentIndex'};
	#
	# Scramble to get the correct version
	if($Pointer < 0)
	{ 
		my $ByteLength = $self->{'Length'};
		$self->DeScramble;
		--$self->{'Rounds'};
		$Pointer = $ByteLength-1;	
		# If the string has a TRIPLE number of bytes, use the
		# last-but-two byte (else previous and next use different 
		# bytes).
		$Pointer -= 2 unless ($ByteLength % 3);
    };
	$self->{'Pointer'} = $Pointer;
	#
	# Result
	return vec(${$self->{'String'}}, $Pointer, 8); 
}

sub string  # (length) -> pointer to random string of length
{
	my $self = shift;
	my $StringLength = shift || 32;
	#
	my $i;
	my $result = "\000" x $StringLength;
        my $PRNGbyte;
	for($i = 0; $i < $StringLength; ++$i)
	{
		#
		# Get the PRNoise as fast as possible!
		# This is inline code for $self->{'NOISE'}->nextByte
		$self->{'Pointer'} += 2;
		++$self->{'CurrentIndex'};
		$PRNGbyte = ($self->{'Pointer'} >= $self->{'Length'})
				? $self->nextByte
				: vec(${$self->{'String'}},
						$self->{'Pointer'}, 8);
		vec($result, $i, 8) = $PRNGbyte;
	};
	# 
	# Return a REFERENCE to the string
	return \$result;
}

# Shuffle a list(-pointer)
sub Shuffle    # (\@list) -> \@list
{
	my $self = shift;
	my $list = shift;
	#
	my $length = scalar(@$list);
	#
	# Shuffle the list
    for($i = ($length-1); $i >= 0; --$i)
    {
       my $j = int($self->fraction * $length);
       my $Value = $$list[$i];
	   $$list[$i] = $$list[$j];
	   $$list[$j] = $Value;
    };
	
	return $list;
}

# Test the integrity of the code
sub Test # () -> 1/0
{
	my $String = "a" x 65;
	my $Char = vec($String, 0, 8);
	my $Test = new KSC_PRNG($String);
	$Test->indexByte(1024);
	$Result1 = $Test->string(32);
	#
	$Test->indexByte(-1);
	$Test->DeRandomize;  # Return to original seed
	$Result2 = $Test->string(32);
	my $i;
	my $n = 0;
	for($i = 0; $i < 32; ++$i)
	{
		return 0 if vec($$Result2, $i, 8) != $Char;
		++$n if vec($$Result1, $i, 8) == $Char;
	};
	return $n < 3 ? 1 : 0;
}

# Read and shuffle all lines from STDIN,
# but only when NOT called inside another program.
unless(caller())
{
	# Initialize a PRNG with the first argument
	my $Seed = shift || "";
	my $PRNG = new KSC_PRNG($Seed);
	
	# Read input
	my @List = <>;
	
	# Shuffle the list
	$PRNG->Shuffle(\@List);
	
	# Print shuffled list
	print STDOUT @List;
};


# Make require happy
1;

