# a script that calls base methylation from bismark_aa --extendedse output
# get the number of 
# at least 10 read coverage 
# at least 20 quality score per C

#use MLDBM;
#use DB_File;
use Getopt::Long;
use List::Util qw(sum);
#use Proc::ProcessTable;

#sub memory_usage {
#  my $t = new Proc::ProcessTable;
#  foreach my $got ( @{$t->table} ) {
#    next if not $got->pid eq $$;
#    return $got->size;
#  }
#}
#my $t = localtime( );warn "$t\n";warn 'memory: '. memory_usage()/1024/1024 ."\n";

# variables
my $minqual=20;
my $mincov =10;
my $phred64= 0;
my $CpGfile ='';
my $CHGfile ='';
my $CHHfile ='';
GetOptions(
	   'read1=s' =>\$read1,
	   'type=s' =>\$type,
	   'nolap!' =>\$nolap,
	   'minqual=i' =>\$minqual,
           'mincov=i' =>\$mincov,
	   'phred64!' =>\$phred64,
	   'CpG=s'=> \$CpGfile,
	   'CHH=s'=> \$CHHfile,
	   'CHG=s'=> \$CHGfile
);
sub usage
{
print STDERR <<EndOfUsage;

Usage: perl $0 [options] input_file >outFile.bed


options:
--read1    : Must be provided at all cases, if given '-' the STDIN will be the input
--type     : one of the following: "single_sam","paired_sam","single_bismark","paired_bismark"
--nolap    : if given and if the input is paired the overlapping paired reads will be ignored
--minqual  : minquality   (default:20)
--mincov   : min coverage (default:10)
--phred64  : quality scores phred64 scale used otherwise phred33 is the default
--CpG      : output filename for CpG methylation scores (if not specified no file is written out)
--CHH      : output filename for CHH methylation scores (if not specified no file is written out)
--CHG      : output filename for CHG methylation scores (if not specified no file is written out)

IMPORTANT:
Files must be sorted based on chr and start of reads. In case of paired-end sam file from bismark, the file still must be sorted in the same way.



EndOfUsage
exit;
} 

if(! defined $type){print "--type argument not supplied\n";usage();}


if(! defined $read1){print "--read1 argument not supplied\n"; usage();}

my $fh;
if($read1 eq "-"){$fh=*STDIN;}
elsif(! -e $read1){
         print "the value of --read1 argument does not point to an existing file\n";
}
else{open ($fh,$read1);}

my %types = map { $_ => 1 }  qw(single_sam paired_sam single_bismark paired_bismark);

if(! exists($types{$type})){
  print "--type argument must be one of the following: 'single_sam','paired_sam','single_bismark','paired_bismark' \n";
  
}

# ARRANGE OFFSET
my $offset=33;
if($phred64){$offset=64;}

# example line from bismark sam output version 0.6.3
# HWI-ST986_0098:1:1101:18264:11272#0/1	0	chr1	497	255	50M	*	0	0	TGGGTTTGATTTGAGGAGAATTGTGTTTCGTTTTTAGAGTATTATCGAAA	CCCFFFFFHHHHHJJJIIJJJJJHHIIJIJHIJJJJJGIDHIJJJIIHJI	NM:i:13	XX:Z:C4C3CC9C4C1C2CC2C6CC1C5	XM:Z:z....x...hx.........x....h.xZ.hh..x......hh.xZ....	XR:Z:CT	XG:Z:CT
#CGGGTCTGACCTGAGGAGAACTGTGCTCCGCCTTCAGAGTACCACCGAAA
#TGGGTTTGATTTGAGGAGAATTGTGTTTCGTTTTTAGAGTATTATCGAAA


#read CATCCTACCAAACTTCGAAATAATCTCCCATATTATAACTTACTACCCCG XX:Z:10GG5G16G11T3	XM:	XR:Z:CT	XG:Z:GA
#gnom CATCCTACCAGGCTTCGGAATAATCTCCCATATTGTAACTTACTACTCCG
#     ..........GG.....G................G...........T...
#     ..........xh....Zx................h..............Z

#C--
# bismark output codes are:
# z unmeth CpG
# Z meth   CpG
# x unmeth CHG
# X meth   CHG
# h unmeth CHH
# H meth   CHH


#open CpG output file and print the header
# print the methylations

if($type eq "single_sam"){
  process_sam($fh,$CpGfile,$CHHfile,$CHGfile,$offset,$mincov,$minqual,0,0);
}
elsif( $type eq "single_bismark"){
  process_single_bismark($fh,$CpGfile,$CHHfile,$CHGfile,$offset,$mincov,$minqual);

}
elsif( $type eq "paired_bismark"){
  die("--paired_bismark option NOT IMPLEMENTED! get a paired sam file and used that as input\n");

}
elsif( $type eq "paired_sam"){
  process_sam($fh,$CpGfile,$CHHfile,$CHGfile,$offset,$mincov,$minqual,1,$nolap);
  #process_single_sam($read2,$CpGfile,$CHHfile,$CHGfile,$offset,0);

}


###  SUBROUTINES ###################

# process a given CG methlation hash
# writes the filter passing CGs to output file
sub processCGmethHash
{
  my($CGmethHash,$out,$mincov)=(@_);

  foreach my $key (keys %{$CGmethHash})
  {
        my($strand,$chr,$loc,$nextBase)=split(/\|/,$key);
	my $noCs=$CGmethHash->{$key}->[0];
	my $noTs=$CGmethHash->{$key}->[1];
	my $noOs=$CGmethHash->{$key}->[2];
	my $Cperc=sprintf("%.2f", 100*$noCs/($noTs+$noCs+$noOs) );
	my $Tperc=sprintf("%.2f", 100*$noTs/($noTs+$noCs+$noOs) );
       	if(($noTs+$noCs)/($noTs+$noCs+$noOs) > 0.9 && ($noTs+$noCs+$noOs)>=$mincov ){
	  print $out join("\t",($chr.".".$loc,$chr,$loc,$strand,$noCs+$noTs+$noOs,$Cperc,$Tperc)  ),"\n"; 
        }
  }
  return 1;
}

# process a given non CG methlation hash
# writes the filter passing Cs to a hash, that hash will be used to calculate conversion rate later on
sub processnonCGmethHash
{
  my($nonCGmethHash,$CTconvArray,$mincov)=(@_);

  foreach my $key (keys %{$nonCGmethHash})
  {
        my($strand,$chr,$loc)=split(/\|/,$key);
	my $noCs=$nonCGmethHash->{$key}->[0];
	my $noTs=$nonCGmethHash->{$key}->[1];
	my $noOs=$nonCGmethHash->{$key}->[2];
	my $Cperc=sprintf("%.2f", 100*$noCs/($noTs+$noCs+$noOs) );
	my $Tperc=sprintf("%.2f", 100*$noTs/($noTs+$noCs+$noOs) );
       	if(($noTs+$noCs)/($noTs+$noCs+$noOs) > 0.95 && ($noTs+$noCs+$noOs)>=$mincov ){
	  #print join("\t",($chr.".".$loc,$chr,$loc,$strand,$noCs+$noTs+$noOs,$Cperc,$Tperc)  ),"\n"; 
	  push @{$CTconvArray->{$strand}},(($noTs*100)/($noTs+$noCs+$noOs));
        }
  }
  return 1;
}



sub processCHmethHash
{
  my($CGmethHash,$out,$mincov)=(@_);

  foreach my $key (keys %{$CGmethHash})
  {
        my($strand,$chr,$loc )=split(/\|/,$key);
	my $noCs=$CGmethHash->{$key}->[0];
	my $noTs=$CGmethHash->{$key}->[1];
	my $noOs=$CGmethHash->{$key}->[2];
	my $Cperc=sprintf("%.2f", 100*$noCs/($noTs+$noCs+$noOs) );
	my $Tperc=sprintf("%.2f", 100*$noTs/($noTs+$noCs+$noOs) );
       	if(($noTs+$noCs)/($noTs+$noCs+$noOs) > 0.9 && ($noTs+$noCs+$noOs)>=$mincov ){
	  print $out join("\t",($chr.".".$loc,$chr,$loc,$strand,$noCs+$noTs+$noOs,$Cperc,$Tperc )  ),"\n"; 
        }
  }
  return 1;
}

# process the methylation call string
sub process_call_string{

  my($mcalls,$i,$key,$CGmethHash,$nonCGmethHash, $CHHmethHash, $CHGmethHash)=@_;

	if( uc($mcalls->[$i]) eq "Z"){ # if genomic base is CpG
	  unless(exists $CGmethHash->{$key}){$CGmethHash->{$key}=[0,0,0];} 
	  if( $mcalls->[$i] eq "Z" ){$CGmethHash->{$key}->[0]++;}          # update Cs
	  elsif( $mcalls->[$i] eq "z"){$CGmethHash->{$key}->[1]++;}        # update Ts
	  else{$CGmethHash->{$key}->[2]++;}                              # update other bases
	}else{                    #if genomic base is non-CpG
	  unless(exists $nonCGmethHash->{$key}){$nonCGmethHash->{$key}=[0,0,0];}

	  my $isCHG;
	  if(uc($mcalls->[$i]) eq "X" && (! exists $CHGmethHash->{$key}) ){
	    $CHGmethHash->{$key}=[0,0,0];
	  }
	  elsif(uc($mcalls->[$i]) eq "H" && (! exists $CHHmethHash->{$key})){
	    $CHHmethHash->{$key}=[0,0,0];
	  }

	  if( $mcalls->[$i] eq "X" )
	  {
	      $nonCGmethHash->{$key}->[0]++;
	      $CHGmethHash->{$key}->[0]++;
	  }
	  elsif($mcalls->[$i] eq "H" )
	  {
	    $nonCGmethHash->{$key}->[0]++;
	    $CHHmethHash->{$key}->[0]++;
	  }
	  elsif( $mcalls->[$i] eq "x" ){
	    $nonCGmethHash->{$key}->[1]++;
	    $CHGmethHash->{$key}->[1]++;

	  }
	  elsif( $mcalls->[$i] eq "h" ){
	    $nonCGmethHash->{$key}->[1]++;
	    $CHHmethHash->{$key}->[1]++;

	  }
	  else{   # this condition will never be used
	    $nonCGmethHash->{$key}->[2]++;
	    if(uc($mcalls->[$i]) eq "X"){	  $CHGmethHash->{$key}->[2]++;}
	    else{	  $CHHmethHash->{$key}->[2]++;}
	  }
	}


}

# get the median value of a given array
# array of numbers
sub median {
    my $rpole = shift;
    my @pole = @$rpole;

    my $ret;

    sort(@pole);

    if( (@pole % 2) == 1 ) {
        $ret = $pole[((@pole+1) / 2)-1];
    } else {
        $ret = ($pole[(@pole / 2)-1] + $pole[@pole / 2]) / 2;
    }

    return $ret;
}

# processed the sam file
sub process_sam{
  my($fh,$CpGfile,$CHHfile,$CHGfile,$offset,$mincov,$minqual,$nolap,$paired)=(@_);

  # check the file status produce flags 
  my $CpGstatus=0;my $CHHstatus=0;my $CHGstatus=0;
  my $out;my $CHHout;my $CHGout;

  if($CpGfile ne ''){
    open ($out,">".$CpGfile);
    print $out join("\t",qw(chrBase chr base strand coverage freqC freqT)),"\n";
    $CpGstatus=1;}
  if($CHHfile ne ''){
    open ($CHHout,">".$CHHfile);
    print $CHHout join("\t",qw(chrBase chr base strand coverage freqC freqT)),"\n";
    $CHHstatus=1;
  }
  if($CHGfile ne '' ){
    open ($CHGout,">".$CHGfile);
    print $CHGout join("\t",qw(chrBase chr base strand coverage freqC freqT)),"\n";
    $CHGstatus=1;
  }

  # check if the file looks like sam
  #read-in file to count C bases
  my %CGmethHash=(); 
  my %nonCGmethHash=(); 
  my %pMeth_nonCG=();
  my %CHHmethHash=();
  my %CHGmethHash=();


  my $lastPos  =-1;
  my $lastChrom="null";
  my $chrPre;
  my $startPre=-1;

  
  while(<$fh>)
  {
    if($_=~/Bismark/){next;} # step over the header line
    if($_=~/^@/){next;} # step over the header line
    # example paired-end reads in SAM format (2 consecutive lines)
    # 1_R1/1	67	5	103172224	255	40M	=	103172417	233	AATATTTTTTTTATTTTAAAATGTGTATTGATTTAAATTT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	NM:i:4	XX:Z:4T1T24TT7	XM:Z:....h.h........................hh.......	XR:Z:CT	XG:Z:CT
    # 1_R1/2	131	5	103172417	255	40M	=	103172224	-233	TATTTTTTTTTAGAGTATTTTTTAATGGTTATTAGATTTT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	NM:i:6	XX:Z:T5T1T9T9T7T3	XM:Z:h.....h.h.........h.........h.......h...	XR:Z:GA	XG:Z:CT
# HWI-ST986_0098:1:1101:18264:11272#0/1	0	chr1	497	255	50M	*	0	0	TGGGTTTGATTTGAGGAGAATTGTGTTTCGTTTTTAGAGTATTATCGAAA	CCCFFFFFHHHHHJJJIIJJJJJHHIIJIJHIJJJJJGIDHIJJJIIHJI	NM:i:13	XX:Z:C4C3CC9C4C1C2CC2C6CC1C5	XM:Z:z....x...hx.........x....h.xZ.hh..x......hh.xZ....	XR:Z:CT	XG:Z:CT

    chomp;
    my @cols   = split(/\t/,$_);
    my $start  = $cols[3];
    my $end    = $start+length($cols[9])-1;
    my $chr    = $cols[2];
    my $methc  = $cols[13]; $methc =~ s/^XM:Z://;
    my @mcalls = split("", $methc); # get the bismark methylation calls
    my @quals  = split("",$cols[10]);# get the quality scores
    my $mrnm   = $cols[6];
    my $mpos   = $cols[7];
    my $isize  = $cols[8];
    my $slen   = length($cols[9]);

    # get strand
    my $strand;
    if($cols[14] eq "XR:Z:CT" && $cols[15] eq "XG:Z:CT" ){$strand="+";}
    elsif($cols[14] eq "XR:Z:CT" && $cols[15] eq "XG:Z:GA"){$strand="-";}
    elsif($cols[14] eq "XR:Z:GA" && $cols[15] eq "XG:Z:CT"){$strand="-";}
    elsif($cols[14] eq "XR:Z:GA" && $cols[15] eq "XG:Z:GA"){$strand="+";}

    # if there is no_overlap trim the mcalls and $quals
    # adjust the start
    if($nolap && ( ($mrnm eq "=") && $paired ) ){
      
      if( ($start+$slen-1)>$mpos){
	if(($mpos-$start)<=0){next;}
	splice @mcalls,($mpos-$start);
	splice @quals,($mpos-$start);
      }

    }



    #checking if the file is sorted
    if( $chr eq $chrPre) {
      if($startPre > $start ){
	die("The sam file is not sorted properly; you can sort the file in unix-like machines using:\n grep -v \'^[[:space:]]*\@\' test.sam | sort -k3,3 -k4,4n  > test.sorted.sam \n");
      }
      $chrPre=$chr;
      $startPre=$start;
    }else{
      $startPre=$start;
      $chrPre=$chr;
    }


    #processes hashes if start-LastPos>100
    if( ($start-$lastPos>100 && $lastPos != -1 ) || ($chr ne $lastChrom && $lastChrom ne "null"  ))
    {
      # if the user wants to write out files write them
      if($CpGstatus){ processCGmethHash(\%CGmethHash,$out,$mincov);}
      if($CHHstatus){ processCHmethHash(\%CHHmethHash,$CHHout,$mincov); }
      if($CHGstatus){ processCHmethHash(\%CHGmethHash,$CHGout,$mincov);}

      processnonCGmethHash(\%nonCGmethHash,\%pMeth_nonCG,$mincov);
      %nonCGmethHash=();
      %CGmethHash=();
      %CHHmethHash=();
      %CHGmethHash=();
    }

    # iterate over the mapped sequence
    for( my $i=0;$i< @quals; $i++) 
    {
	if ( (ord($quals[$i])-$offset) < $minqual || ( $mcalls[$i] eq ".") ){next;}
	#if(( $gbases[$i] eq "C" && $i==(scalar @quals)) && $gbases[$i-1].$gbases[$i].$gbases[$i+1].$gbases[$i+2] eq "CCGG" ){next;} #if last base is a C and it is a part of CCGG motif, don't call for meth
	my $key;# initilaize the hash key
	if($strand eq "+"){$key=join("|",("F",$chr,$start+$i));}else{$key=join("|",("R",$chr,$start+$i));}

	process_call_string(\@mcalls,$i,$key, \%CGmethHash, \%nonCGmethHash,  \%CHHmethHash,  \%CHGmethHash);


    }
    $lastPos=$end;
    $lastChrom=$chr;
  }
  close $fh;

  if($CpGstatus){ processCGmethHash(\%CGmethHash,$out,$mincov);}
  if($CHHstatus){ processCHmethHash(\%CHHmethHash,$CHHout,$mincov); }
  if($CHGstatus){ processCHmethHash(\%CHGmethHash,$CHGout,$mincov);}
  processnonCGmethHash(\%nonCGmethHash,\%pMeth_nonCG,$mincov);

  #close $out;


  # get the conversion rate and write it out!!

  my $numF=scalar @{$pMeth_nonCG{"F"}};
  my $numR=scalar @{$pMeth_nonCG{"R"}};

  if($numF==0 && $numR==0){
    #if($CpGstatus){unlink($CpGfile);}
    #if($CHHstatus){unlink($CHHfile);}
    #if($CHGstatus){unlink($CHGfile);}
    die("\nnot enough alignments that pass coverage and phred score thresholds to calculate conversion rates\n EXITING....\n\n");}

  my $AvFconvRate=0;
  my $AvRconvRat=0;
  my $medFconvRate=0;
  my $medRconvRate=0;
  if($numF>0){ $AvFconvRate=(sum(@{$pMeth_nonCG{"F"}}))/$numF; }
 if($numR>0){ $AvRconvRate=(sum(@{$pMeth_nonCG{"R"}}))/$numR; }
  my $AvconvRate =(sum(@{$pMeth_nonCG{"F"}})+ sum( @{$pMeth_nonCG{"R"}}))/($numF+$numR);

  my @allesSchon;push @allesSchon,@{$pMeth_nonCG{"F"}},@{$pMeth_nonCG{"R"}};
  if($numF>0){ $medFconvRate=median($pMeth_nonCG{"F"}); }
  if($numR>0){ $medRconvRate=median($pMeth_nonCG{"R"}); }
  my $medconvRate =median(\@allesSchon);

  my $totCpG=(scalar(@allesSchon));

  my $res="";
  $res .= "total otherC considered (>95% C+T): $totCpG\n";
  $res .= "average conversion rate = $AvconvRate\n"; 
  $res .= "median conversion rate = $medconvRate\n\n"; 

  $res .= "total otherC considered (Forward) (>95% C+T): $numF\n";
  $res .= "average conversion rate (Forward) = $AvFconvRate\n"; 
  $res .= "median conversion rate (Forward) = $medFconvRate\n\n"; 

  $res .= "total otherC considered (Reverse) (>95% C+T): $numR\n";
  $res .= "average conversion rate (Reverse) = $AvRconvRate\n"; 
  $res .= "median conversion rate (Reverse) = $medRconvRate\n"; 

  #open (my $hd,">".$prefix."_conversionRate.txt");
  #print $hd $res;
  #close $hd;

  print $res;


}





# processed the sam file
sub process_single_bismark{
  my($fh,$CpGfile,$CHHfile,$CHGfile,$offset,$mincov,$minqual)=(@_);

  # check the file status produce flags 
  my $CpGstatus=0;
  my $CHHstatus=0;
  my $CHGstatus=0;
  my $out;
  my $CHHout;
  my $CHGout;

  if($CpGfile ne ''){

    open ($out,">".$CpGfile);
    print $out join("\t",qw(chrBase chr base strand coverage freqC freqT)),"\n";
    $CpGstatus=1;
  }

  if($CHHfile ne ''){
    open ($CHHout,">".$CHHfile);
    print $CHHout join("\t",qw(chrBase chr base strand coverage freqC freqT)),"\n";
    $CHHstatus=1;
  }

  if($CHGfile ne '' ){
    open ($CHGout,">".$CHGfile);
    print $CHGout join("\t",qw(chrBase chr base strand coverage freqC freqT)),"\n";
    $CHGstatus=1;
  }

  # check if the file looks like sam
  #read-in file to count C bases
  my %CGmethHash=(); 
  my %nonCGmethHash=(); 
  my %pMeth_nonCG=();
  my %CHHmethHash=();
  my %CHGmethHash=();




  my $chrPre;
  my $startPre=-1;
  my $lastPos  =-1;
  my $lastChrom="null";
  while(<$fh>)
  {
    if($_=~/Bismark/){next;} # step over the header line
    if($_=~/^@/){next;} # step over the header line

    chomp;
    my @cols   = split(/\t/,$_);
    my $start  = $cols[3];my $end   =$cols[4];my $strand=$cols[1];my $chr   =$cols[2];
    my @mcalls = split("",$cols[7]); # get the bismark methylation calls
    my @gbases = split("",$cols[6]);# get the genomic bases
    my @quals  = split("",$cols[10]);# get the quality scores


    #checking if the file is sorted
    if( $chr eq $chrPre) {
      if($startPre > $start ){
	die("The sam file is not sorted properly; you can sort the file in unix-like machines using:\n grep -v \'^[[:space:]]*\@\' test.sam | sort -k3,3 -k4,4n  > test.sorted.sam \n");
      }
      $chrPre=$chr;
      $startPre=$start;
    }else{
      $startPre=$start;
      $chrPre=$chr;
    }


    #processes hashes if start-LastPos>100
    if( ($start-$lastPos>100 && $lastPos != -1 ) || ($chr ne $lastChrom && $lastChrom ne "null"  ))
    {
      # if the user wants to write out files write them
      if($CpGstatus){ processCGmethHash(\%CGmethHash,$out,$mincov);}
      if($CHHstatus){ processCHmethHash(\%CHHmethHash,$CHHout,$mincov); }
      if($CHGstatus){ processCHmethHash(\%CHGmethHash,$CHGout,$mincov);}

      processnonCGmethHash(\%nonCGmethHash,\%pMeth_nonCG,$mincov);
      %nonCGmethHash=();
      %CGmethHash=();
      %CHHmethHash=();
      %CHGmethHash=();
    }

      # iterate over the mapped sequence
      for( my $i=0;$i< @quals; $i++) 
      {
	if ( (ord($quals[$i])-$offset) < $minqual || ( $mcalls[$i] eq ".") ){next;}
	if(( $gbases[$i] eq "C" && $i==(scalar @quals)) && $gbases[$i-1].$gbases[$i].$gbases[$i+1].$gbases[$i+2] eq "CCGG" ){next;} #if last base is a C and it is a part of CCGG motif, don't call for meth
	my $key;# initilaize the hash key
	if($strand eq "+"){$key=join("|",("F",$chr,$start+$i));}else{$key=join("|",("R",$chr,$end-$i));}
	process_call_string(\@mcalls,$i,$key, \%CGmethHash, \%nonCGmethHash,  \%CHHmethHash,  \%CHGmethHash);

      }
      $lastPos=$end;
      $lastChrom=$chr;
  }
  close $fh;
  close $out;


  # get the conversion rate and write it out!!

  my $numF=scalar @{$pMeth_nonCG{"F"}};
  my $numR=scalar @{$pMeth_nonCG{"R"}};

  if($numF==0 || $numR==0){
    #if($CpGstatus){unlink($CpGfile);}
    #if($CHHstatus){unlink($CHHfile);}
    #if($CHGstatus){unlink($CHGfile);}
    die("\nnot enough aligments that pass coverage and phred score thresholds to calculate conversion rates\nEXITING....\n\n");}

  my $AvFconvRate=(sum(@{$pMeth_nonCG{"F"}}))/$numF;
  my $AvRconvRate=(sum(@{$pMeth_nonCG{"R"}}))/$numR;
  my $AvconvRate =(sum(@{$pMeth_nonCG{"F"}})+ sum( @{$pMeth_nonCG{"R"}}))/($numF+$numR);

  my @allesSchon;push @allesSchon,@{$pMeth_nonCG{"F"}},@{$pMeth_nonCG{"R"}};
  my $medFconvRate=median($pMeth_nonCG{"F"});
  my $medRconvRate=median($pMeth_nonCG{"R"});
  my $medconvRate =median(\@allesSchon);

  my $totCpG=(scalar(@allesSchon));

  my $res="";
  $res .= "total otherC considered (>95% C+T): $totCpG\n";
  $res .= "average conversion rate = $AvconvRate\n"; 
  $res .= "median conversion rate = $medconvRate\n\n"; 

  $res .= "total otherC considered (Forward) (>95% C+T): $numF\n";
  $res .= "average conversion rate (Forward) = $AvFconvRate\n"; 
  $res .= "median conversion rate (Forward) = $medFconvRate\n\n"; 

  $res .= "total otherC considered (Reverse) (>95% C+T): $numR\n";
  $res .= "average conversion rate (Reverse) = $AvRconvRate\n"; 
  $res .= "median conversion rate (Reverse) = $medRconvRate\n"; 

  #open (my $hd,">".$prefix."_conversionRate.txt");
  #print $hd $res;
  #close $hd;

  print $res;


}


sub process_paired_bismark{
  my($input_file,$CpGfile,$CHHfile,$CHGfile,$offset,$mincov,$minqual,$ignore)=(@_);
}



