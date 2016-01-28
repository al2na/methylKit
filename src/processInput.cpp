#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <deque>
#include <string>
#include <map>
#include <regex>
//#include <Rcpp.h>

bool String2Int(const std::string& str, int& result)
{
    try
    {
        std::size_t lastChar;
        result = std::stoi(str, &lastChar, 10);
        return lastChar == str.size();
    }
    catch (std::invalid_argument&)
    {
        return false;
    }
    catch (std::out_of_range&)
    {
        return false;
    }
}


//split function to seperate strings by delimiter 
// taken from http://stackoverflow.com/a/236803
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}



//###  SUBROUTINES ###################

//# process a given CG methlation hash
//# writes the filter passing CGs to output file
int processCGmethHash(std::map<std::string,std::vector<int> > &CGmethHash,std::FILE* out, int &mincov) {

  // iterate over the map/hash
  for ( std::map<std::string,std::vector<int> >::iterator iter=CGmethHash.begin(); iter!=CGmethHash.end(); ++iter ) {
    // get keys and values pair
    std::string key = iter->first; 
    std::vector<int> value =  iter->second;
    // split key to get strand, chr, and locus
    std::vector<std::string> temp = split(key,'|');
    std::string strand = temp[0] , chr =  temp[1], loc = temp[2];
    //int loc = std::atoi(temp[2]); //nextBase; 
    
    int noCs= value[0];
    int noTs= value[1];
    int noOs= value[2];
    double Cperc= 100*noCs/(double)(noTs + noCs + noOs) ;
    double Tperc= 100*noTs/(double)(noTs + noCs + noOs) ;
    if((( noTs + noCs)/(double)(noTs + noCs + noOs) > 0.9) && ((noTs+noCs+noOs)>=mincov) ){
      // print to file : "chr.loc   chr   loc   strand   totalbases   %C   %T"
      fprintf(out, "%s.%s\t%s\t%s\t%s\t%i\t%.2f\t%.2f\n", 
              chr.c_str(), loc.c_str(), chr.c_str(), loc.c_str(), strand.c_str(), noCs+noTs+noOs, Cperc,Tperc); 
    }
  }
return 0;
}





//# process a given non CG methlation hash
//# writes the filter passing Cs to a hash of arrays, that hash will be used to calculate conversion rate later on
//# one might use binning for storing hash arrays to decrease the memory load for median calcuation


int processnonCGmethHash ( std::map<std::string, std::vector<int> > &nonCGmethHash,std::map<std::string,std::map<std::string, double> >  &CTconvArray, int &mincov) { 

  // iterate over map/hash
  for ( std::map<std::string,std::vector<int> >::iterator iter=nonCGmethHash.begin(); iter!=nonCGmethHash.end(); ++iter ) {
    // get keys and values pair
    std::string key = iter->first; 
    std::vector<int> value =  iter->second;
    // split key to get strand, chr, and locus
    std::vector<std::string> temp = split(key,'|');
    std::string strand = temp[0] , chr =  temp[1], loc = temp[2];
    //int loc = std::atoi(temp[2]); //nextBase; 
    
    int noCs= value[0];
    int noTs= value[1];
    int noOs= value[2];
    // double Cperc= 100*noCs/(double)(noTs + noCs + noOs) ;
    // double Tperc= 100*noTs/(double)(noTs + noCs + noOs) ;
    if((( noTs + noCs)/(double)(noTs + noCs + noOs) > 0.9) && ((noTs+noCs+noOs)>=mincov) ){
      //// print to file : "chr.loc   chr   loc   strand   totalbases   %C   %T"
      //fprintf(out, "%s.%s\t%s\t%s\t%s\t%i\t%.2f\t%.2f\n", 
              //chr.c_str(), loc.c_str(), chr.c_str(), loc.c_str(), strand.c_str(), noCs+noTs+noOs, Cperc,Tperc); 
    //}
    //   push @{$CTconvArray->{$strand}},(($noTs*100)/($noTs+$noCs+$noOs));
      //CTconvArray[strand]["median"] += ((noTs*100)/(noTs+noCs+noOs));
    }
    }
    return 1;
    }

/** process a given non CG methlation hash
/ writes the filter passing Cs to a hash, that hash will be used to calculate conversion rate later on
/ hash only contains total number of non-CpGs Cs and total methylation value to calculate average conversion rate
/ using this we won't be able to calculate median conversion rate, just the mean conversion rate
 this is better for memory management **/

int processnonCGmethHash2 ( std::map<std::string, std::vector<int> > &nonCGmethHash,std::map<std::string,std::map<std::string, double> >  &CTconvArray, int &mincov) { 

  // iterate over map/hash
  for ( std::map<std::string,std::vector<int> >::iterator iter=nonCGmethHash.begin(); iter!=nonCGmethHash.end(); ++iter ) {
    // get keys and values pair
    std::string key = iter->first; 
    std::vector<int> value =  iter->second;
    // split key to get strand, chr, and locus
    std::vector<std::string> temp = split(key,'|');
    std::string strand = temp[0] , chr =  temp[1], loc = temp[2];
    //int loc = std::atoi(temp[2]); //nextBase; 
    
    int noCs= value[0];
    int noTs= value[1];
    int noOs= value[2];
    // double Cperc= 100*noCs/(double)(noTs + noCs + noOs) ;
    // double Tperc= 100*noTs/(double)(noTs + noCs + noOs) ;
    if((( noTs + noCs)/(double)(noTs + noCs + noOs) > 0.9) && ((noTs+noCs+noOs)>=mincov) ){
    //print join("\t",($chr.".".$loc,$chr,$loc,$strand,$noCs+$noTs+$noOs,$Cperc,$Tperc)  ),"\n"; 
    //push @{$CTconvArray->{$strand}},(($noTs*100)/($noTs+$noCs+$noOs));
    CTconvArray[strand]["num"]++;
    CTconvArray[strand]["total"] += ((noTs*100)/(noTs+noCs+noOs));
}
}
return 0;
}


int processCHmethHash ( std::map<std::string, std::vector<int> > &CGmethHash, std::FILE *out, int &mincov) { 

  // iterate over the map/hash
  for ( std::map<std::string,std::vector<int> >::iterator iter=CGmethHash.begin(); iter!=CGmethHash.end(); ++iter ) {
    // get keys and values pair
    std::string key = iter->first; 
    std::vector<int> value =  iter->second;
    // split key to get strand, chr, and locus
    std::vector<std::string> temp = split(key,'|');
    std::string strand = temp[0] , chr =  temp[1], loc = temp[2];
    //int loc = std::atoi(temp[2]); //nextBase; 
    
    int noCs= value[0];
    int noTs= value[1];
    int noOs= value[2];
    double Cperc=  100*noCs/(double) (noTs + noCs + noOs) ;
    double Tperc=  100*noTs/(double) (noTs + noCs + noOs) ;
    if((( noTs + noCs)/(double)(noTs + noCs + noOs) > 0.9) && ((noTs+noCs+noOs)>=mincov) ){
      // print to file : "chr.loc   chr   loc   strand   totalbases   %C   %T"
      fprintf(out, "%s.%s\t%s\t%s\t%s\t%i\t%.2f\t%.2f\n", 
              chr.c_str(), loc.c_str(), chr.c_str(), loc.c_str(), strand.c_str(), noCs+noTs+noOs, Cperc,Tperc); 
    }
  }
return 0;
}

//# process the methylation call string
void process_call_string (std::string &mcalls, int &i, std::string &key, std::map<std::string, std::vector<int> > &CGmethHash, 
                          std::map<std::string, std::vector<int> > &nonCGmethHash, std::map<std::string, std::vector<int> > &CHHmethHash, 
                          std::map<std::string, std::vector<int> > &CHGmethHash) {

  std::pair<std::map<std::string, std::vector<int> >::iterator,bool> ret;


  if( std::toupper( mcalls[i]) == 'Z'){ // if genomic base is CpG
    // insert new key into map if not existent, else return iterator to existing key and change nothing
    CGmethHash.insert(std::pair<std::string, std::vector<int> >(key, std::vector<int>(3,0)));

    if( mcalls[i] == 'Z' )       
      { CGmethHash[key][0]++; }       // update Cs
    else if( mcalls[i] ==  'z') 
      { CGmethHash[key][1]++; }        // update Ts
    else                                    
      { CGmethHash[key][2]++; }        // update other bases
  } else{   //if genomic base is non-CpG
    // insert new key into map if not existent, else return iterator to existing key and change nothing
    nonCGmethHash.insert(std::pair<std::string, std::vector<int> >(key, std::vector<int>(3,0)));
    if( std::toupper( mcalls[i]) == 'X' ) {
        // insert new key into map if not existent, else return iterator to existing key and change nothing
        CHGmethHash.insert(std::pair<std::string, std::vector<int> >(key, std::vector<int>(3,0)));
    }
    else if( std::toupper( mcalls[i] ) == 'H') {
        // insert new key into map if not existent, else return iterator to existing key and change nothing
        CHHmethHash.insert(std::pair<std::string, std::vector<int> >(key, std::vector<int>(3,0)));
    }
    
    if( mcalls[i] == 'X'){
      nonCGmethHash[key][0]++;
      CHGmethHash[key][0]++;
    }
    else if( mcalls[i] == 'H') {
      nonCGmethHash[key][0]++;
      CHHmethHash[key][0]++;
    }
    else if( mcalls[i] == 'x') {
      nonCGmethHash[key][1]++;
      CHGmethHash[key][1]++;
    }
    else if( mcalls[i] == 'h') {
      nonCGmethHash[key][1]++;
      CHHmethHash[key][1]++;
      
    }
    // else {   // this condition will never be used
    //   nonCGmethHash[key][2]++;
    //   if( std::toupper(*mcalls[i]) =='X') { CHGmethHash[key][2]++; }
    //   else { CHHmethHash[key][2]++;}
    // }
  }
}

// get the median value of a given array
// array of numbers
double median(std::vector<double> vec) {
  
  std::nth_element(vec.begin(), vec.begin() + vec.size()/2, vec.end());
  
  if ( vec.size() % 2 == 0)
    return (vec[vec.size()/2 - 1] + vec[vec.size()/2]) /2.0;
  else 
    return vec[vec.size()/2];
  
  
}


// processes the cigar string and remove and delete elements from mcalls and quality scores
void processCigar ( std::string cigar, std::string &methc, std::string &qual) {
  int position = 0, len;
  std::smatch m;
  std::string cigar_part, insertion;


  std::deque<int> insPos; // location of the insertions
  std::deque<int> insLen; // location of the insert lengths
  while (!cigar.empty()){
    if(std::regex_search(cigar, m, std::regex ("^[0-9]+[MIDS]"))) {
      cigar_part = m.str();
      if (std::regex_search(cigar_part, m, std::regex ("(\\d+)M"))) { // count matches
        position += std::stoi(m[1].str());
      } 
      else if (std::regex_search(cigar_part, m, std::regex ("(\\d+)I"))) { // count inserts
        len = std::stoi(m[1].str());
        insertion = std::string ( len ,'-');
        insPos.push_front(position); 
        insLen.push_front(len); 
        
        position += len;
      } 
      else if (std::regex_search(cigar_part, m, std::regex ("(\\d+)D"))) { // count deletions
        len = std::stoi(m[1].str());
        insertion = std::string(len, '.');
        methc.insert(position, insertion);
        qual.insert(position, insertion);
        
        position += len;
      } 
      else if (std::regex_search(cigar_part, std::regex ("\\d+)S"))) { 
        //#############################
        // die "Not ready for this!\n";   
        // ###########################

      }
      cigar.erase(0,cigar_part.length());
    } else {
        // #############################
        // die "Unexpected cigar: $id $cigar\n";   
        // ###########################
    }
  }
  
  // if there are insertions remove it from the mcalls and 
  if(insPos.size()>0){
    for(int i=0;i < insPos.size();i++){
      methc.erase(insPos[i],insLen[i]);
      qual.erase(insPos[i],insLen[i]);
    }
  }
  
}

// processed the sam file
int process_sam ( std::istream &fh, char* CpGfile, char* CHHfile, char* CHGfile, int &offset, int &mincov, int &minqual, int &nolap, int &paired) {

  
// check the file status produce flags 
  int CpGstatus = 0, CHHstatus = 0, CHGstatus = 0;
  FILE *out, *CHHout, *CHGout;
  
  
  if( (CpGfile!=NULL) &&  (std::strlen(CpGfile) != 0)) {
    out = fopen(CpGfile,"w");
    fprintf (out, "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n");
    CpGstatus=1;
  }
  if((CHHfile!=NULL) &&  (std::strlen(CHHfile) != 0)) {
    CHHout = fopen(CHHfile,"w");
    fprintf (CHHout, "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n");
    CHHstatus=1;
  }
  if((CHGfile!=NULL) &&  (std::strlen(CHGfile) != 0)) {
    CHGout = fopen(CHGfile,"w");
    fprintf (CHGout, "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n");
    CHGstatus=1;
  }
  
  //# check if the file looks like sam
  //#read-in file to count C bases
  std::map<std::string, std::vector<int> > CGmethHash; 
  std::map<std::string, std::vector<int> > nonCGmethHash; 
  std::map<std::string, std::vector<int> > CHHmethHash;
  std::map<std::string, std::vector<int> > CHGmethHash;
  std::map<std::string,std::map<std::string, double> > pMeth_nonCG;

  
  int lastPos  =-1, startPre = -1 ;
  std::string lastChrom="null", chrPre;

  std::string line;  
  
  while(std::getline(fh, line))
  {
    
    //std::cout << line << std::endl;
    if(std::regex_search(line, std::regex ("Bismark")))  {std::getline(fh, line);}  // step over the header line
    if(std::regex_search(line, std::regex ("^@")))       {std::getline(fh, line);} // step over the header line
    /** example paired-end reads in SAM format (2 consecutive lines)
    # 1_R1/1	67	5	103172224	255	40M	=	103172417	233	AATATTTTTTTTATTTTAAAATGTGTATTGATTTAAATTT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	NM:i:4	XX:Z:4T1T24TT7	XM:Z:....h.h........................hh.......	XR:Z:CT	XG:Z:CT
    # 1_R1/2	131	5	103172417	255	40M	=	103172224	-233	TATTTTTTTTTAGAGTATTTTTTAATGGTTATTAGATTTT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	NM:i:6	XX:Z:T5T1T9T9T7T3	XM:Z:h.....h.h.........h.........h.......h...	XR:Z:GA	XG:Z:CT
    # HWI-ST986_0098:1:1101:18264:11272#0/1	0	chr1	497	255	50M	*	0	0	TGGGTTTGATTTGAGGAGAATTGTGTTTCGTTTTTAGAGTATTATCGAAA	CCCFFFFFHHHHHJJJIIJJJJJHHIIJIJHIJJJJJGIDHIJJJIIHJI	NM:i:13	XX:Z:C4C3CC9C4C1C2CC2C6CC1C5	XM:Z:z....x...hx.........x....h.xZ.hh..x......hh.xZ....	XR:Z:CT	XG:Z:CT
    **/




    std::vector<std::string> cols = split(line,'\t');
    int start;                     
    if(!String2Int(cols[3],start)) { std::cout << "Error from String2Int" << std::endl;return -1;}
    int end                       = start + cols[9].length()-1;
    std::string chr               = cols[2];
    std::string cigar             = cols[5];
    std::string methc             = cols[13]; 
    methc.erase(methc.begin(),methc.begin()+ 5 ); //  remove "XM:Z:"
    std::string qual              = cols[10];
    std::string mrnm              = cols[6];
    int mpos;                     // = std::stoi(cols[7]);
    int isize;                    // = std::stoi(cols[8]);
    if(!String2Int(cols[7],mpos)) { std::cout << "Error from String2Int" << std::endl;return -1;}
    if(!String2Int(cols[8],isize)) { std::cout << "Error from String2Int" << std::endl;return -1;}
    
    
    // process cigar string to get indels
    if( std::regex_search(cigar, std::regex ("[DI]"))) {
      processCigar( cigar, methc, qual);
    }
    std::string mcalls = methc; // get the bismark methylation calls
    std::string quals  = qual;  // get the quality scores
    int slen   = mcalls.length(); // aligned sequence length
    
    // get strand
    char strand;
    if( (cols[14] == "XR:Z:CT") && (cols[15] == "XG:Z:CT") ) {strand='+';} // original top strand
    else if( (cols[14] == "XR:Z:CT") && (cols[15] == "XG:Z:GA") ) {strand='-';} // original bottom strand
    else if( (cols[14] == "XR:Z:GA") && (cols[15] == "XG:Z:CT") ) {strand='+';} // complementary to original top strand, bismark says - strand to this
    else if( (cols[14] == "XR:Z:GA") && (cols[15] == "XG:Z:GA") ) {strand='-';} // complementary to original bottom strand, bismark says + strand to this
    
    // if there is no_overlap trim the mcalls and $quals
    // adjust the start
    if( nolap && ( ( mrnm == "=") && paired ) ){
      if( ( start + slen - 1) > mpos) {
        if( ( mpos - start ) > 0 ) { //{continue;}
          mcalls.erase(mpos-start ,std::string::npos);
          quals.erase(mpos-start, std::string::npos);
        }
      }
    }
    
    
    
    //checking if the file is sorted
    if( chr == chrPre) {
      if( startPre > start ) {
        // ####################### 
        //die("The sam file is not sorted properly; you can sort the file in unix-like machines using:\n grep -v \'^[[:space:]]*\@\' test.sam | sort -k3,3 -k4,4n  > test.sorted.sam \n");
        // ########################
      }
      chrPre=chr;
      startPre=start;
    } else {
      startPre=start;
      chrPre=chr;
    }
    
    
    //processes hashes if start-LastPos>100
    if( ( (start- lastPos > 100) && (lastPos != -1 )) || ( (chr != lastChrom) && (lastChrom != "null")  ))
    {
      // if the user wants to write out files write them
      if(CpGstatus){ processCGmethHash(CGmethHash,out,mincov);}
      if(CHHstatus){ processCHmethHash(CHHmethHash,CHHout,mincov); }
      if(CHGstatus){ processCHmethHash(CHGmethHash,CHGout,mincov);}
      
      //processnonCGmethHash(\%nonCGmethHash,\%pMeth_nonCG,$mincov);
      processnonCGmethHash2(nonCGmethHash,pMeth_nonCG,mincov);
      
      nonCGmethHash.clear();
      CGmethHash.clear();
      CHHmethHash.clear();
      CHGmethHash.clear();
    }
    
    // iterate over the mapped sequence
    for( int i=0; i < quals.length(); i++) 
    {
      if ( ( ( (int) quals[i] - offset ) <  minqual ) || ( mcalls[i] == '.') ){ continue;}
      std::string key; // initialize the hash key
      if( strand == '+') { key = "F|"+ chr+"|"+std::to_string(start+i); }
      else { key = "R|"+ chr+"|"+std::to_string(start+i); }

      process_call_string(mcalls,i,key, CGmethHash, nonCGmethHash, CHHmethHash, CHGmethHash);
      
      
    }
    lastPos=end;
    lastChrom=chr;
  }
  //fh.close();
  
  if(CpGstatus) { processCGmethHash(CGmethHash,out,mincov);     std::fclose(out); }
  if(CHHstatus) { processCHmethHash(CHHmethHash,CHHout,mincov); std::fclose(CHHout); }
  if(CHGstatus) { processCHmethHash(CHGmethHash,CHGout,mincov); std::fclose(CHGout); }
  //processnonCGmethHash(nonCGmethHash,pMeth_nonCG,mincov);
  processnonCGmethHash2(nonCGmethHash,pMeth_nonCG,mincov);


  // get the conversion rate and write it out!!
  
  
  int numF= pMeth_nonCG["F"]["num"];
  int numR= pMeth_nonCG["R"]["num"];
  //std::printf( "%i %i\n", numF, numR);
  // if( (numF == 0) && (numR == 0)) {
  //   if(CpGstatus){std::remove(CpGfile);}
  //   if(CHHstatus){std::remove(CHHfile);}
  //   if(CHGstatus){std::remove(CHGfile);}
  //   // #################
  //   //die("\nnot enough alignments that pass coverage and phred score thresholds to calculate conversion rates\n EXITING....\n\n");
  //   // ################
  // }


  double AvFconvRate  = 0 , AvRconvRate  = 0;
  //int medFconvRate = 0,  medRconvRate = 0;
  if(numF > 0) { AvFconvRate = pMeth_nonCG["F"]["total"]/numF; }
  if(numR > 0) { AvRconvRate = pMeth_nonCG["R"]["total"]/numR; }
  double AvconvRate =(pMeth_nonCG["F"]["total"] + pMeth_nonCG["R"]["total"])/(numF+numR);
  
  //std::vector<double> allesSchon; allesSchon.push_back(pMeth_nonCG["F")); allesSchon.push_back(pMeth_nonCG["R"));
  //if( numF>0 ) { medFconvRate=median(pMeth_nonCG["F"]); }
  //if( numR>0 ) { medRconvRate=median(pMeth_nonCG["R"]); }
  //double medconvRate = median(\@allesSchon);
  
  //int totCs=allesSchon.size();
  int totCs = numF + numR;

   
  std::cout 
    <<  "total otherC considered (>95%% C+T): "           <<   std::setprecision(9) << totCs       << "\n"
    <<  "average conversion rate = "                      <<   std::setprecision(9)<< AvconvRate  << "\n"
      //"median conversion rate = %.2f\n\n";              //$medconvRate
    
    <<  "total otherC considered (Forward) (>95%% C+T): " <<   std::setprecision(9)<< numF        << "\n"
    <<  "average conversion rate (Forward) = "            <<   std::setprecision(9)<< AvFconvRate << "\n"
      //"median conversion rate (Forward) = %.2f\n\n"        //$medFconvRate 
    
    <<  "total otherC considered (Reverse) (>95%% C+T): " <<   std::setprecision(9)<< numR        << "\n"
    <<  "average conversion rate (Reverse) = "            <<   std::setprecision(9)<< AvRconvRate << "\n";
  
  //open (my $hd,">".$prefix."_conversionRate.txt");
  //print $hd $res;
  //close $hd;
  
  return 0;
  
}





// processed the sam file
int process_single_bismark (std::istream& fh, char* CpGfile, char* CHHfile, char* CHGfile, int &offset, int &mincov, int &minqual ) {

  
// check the file status produce flags 
  int CpGstatus = 0, CHHstatus = 0, CHGstatus = 0;
  FILE *out, *CHHout, *CHGout;
  
  
  if( (CpGfile!=NULL) &&  (std::strlen(CpGfile) != 0)) {
    out = fopen(CpGfile,"w");
    fprintf (out, "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n");
    CpGstatus=1;
  }
  if((CHHfile!=NULL) &&  (std::strlen(CHHfile) != 0)) {
    CHHout = fopen(CHHfile,"w");
    fprintf (CHHout, "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n");
    CHHstatus=1;
  }
  if((CHGfile!=NULL) &&  (std::strlen(CHGfile) != 0)) {
    CHGout = fopen(CHGfile,"w");
    fprintf (CHGout, "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n");
    CHGstatus=1;
  }
  
  //# check if the file looks like sam
  //#read-in file to count C bases
  std::map<std::string, std::vector<int> > CGmethHash; 
  std::map<std::string, std::vector<int> > nonCGmethHash; 
  std::map<std::string, std::vector<int> > CHHmethHash;
  std::map<std::string, std::vector<int> > CHGmethHash;
  std::map<std::string,std::map<std::string, double> > pMeth_nonCG;
  
  
  
  
  int lastPos  =-1, startPre = -1 ;
  std::string lastChrom="null", chrPre;

  std::string line;  
  
  while(fh)
  {
    std::getline(fh, line);
    if(std::regex_search(line, std::regex ("Bismark")))  {std::getline(fh, line);}  // step over the header line
    if(std::regex_search(line, std::regex ("^@")))       {std::getline(fh, line);} // step over the header line
    
    std::vector<std::string> cols = split(line,'\t');
    int start , end;                    
    if(!String2Int(cols[3],start)) { std::cout << "Error from String2Int" << std::endl;return -1;}
    if(!String2Int(cols[4],end)) { std::cout << "Error from String2Int" << std::endl;return -1;}
    char strand                   = cols[1][0];
    std::string chr               = cols[2]; 
    std::string qual              = cols[10];
    std::string mcalls            = cols[7];  // get the bismark methylation calls
    std::string gbases            = cols[6];  // get the genomic bases
    std::string quals             = cols[10]; // get the quality scores

    
    
    //checking if the file is sorted
    if( chr == chrPre) {
      if( startPre > start ) {
        // ####################### 
        //die("The sam file is not sorted properly; you can sort the file in unix-like machines using:\n grep -v \'^[[:space:]]*\@\' test.sam | sort -k3,3 -k4,4n  > test.sorted.sam \n");
        // ########################
      }
      chrPre=chr;
      startPre=start;
    } else {
      startPre=start;
      chrPre=chr;
    }
    
    
    //processes hashes if start-LastPos>100
    if( ( (start- lastPos > 100) && (lastPos != -1 )) || ( (chr != lastChrom) && (lastChrom != "null")  ))
    {
      // if the user wants to write out files write them
      if(CpGstatus){ processCGmethHash(CGmethHash,out,mincov);}
      if(CHHstatus){ processCHmethHash(CHHmethHash,CHHout,mincov); }
      if(CHGstatus){ processCHmethHash(CHGmethHash,CHGout,mincov);}
      
      //processnonCGmethHash(\%nonCGmethHash,\%pMeth_nonCG,$mincov);
      processnonCGmethHash2(nonCGmethHash,pMeth_nonCG,mincov);
      
      nonCGmethHash.clear();
      CGmethHash.clear();
      CHHmethHash.clear();
      CHGmethHash.clear();
    }
    
    // iterate over the mapped sequence
    for( int i=0; i < quals.length(); i++) 
    {
      if ( ( ( (int) quals[i] - offset ) <  minqual ) || ( mcalls[i] == '.') ){ continue;}
      //if last base is a C and it is a part of CCGG motif, don't call for meth
      if( ( (gbases[i] == 'C') && (i==quals.length()) ) && ( gbases.substr(i-1,4) == "CCGG" ) ) { continue;} 
      std::string key; // initilaize the hash key
      if( strand == '+') { key = "F|"+ chr+"|"+std::to_string(start+i); }
      else { key = "R|"+ chr+"|"+std::to_string(start+i); }
      
      process_call_string(mcalls, i,key, CGmethHash, nonCGmethHash, CHHmethHash, CHGmethHash);
      
      
    }
    lastPos=end;
    lastChrom=chr;
  }
  //close $fh;
  fclose(out);


  if(CpGstatus) { processCGmethHash(CGmethHash,out,mincov);     std::fclose(out); }
  if(CHHstatus) { processCHmethHash(CHHmethHash,CHHout,mincov); std::fclose(CHHout); }
  if(CHGstatus) { processCHmethHash(CHGmethHash,CHGout,mincov); std::fclose(CHGout); }
  //processnonCGmethHash(nonCGmethHash,pMeth_nonCG,mincov);
  processnonCGmethHash2(nonCGmethHash,pMeth_nonCG,mincov);


  // get the conversion rate and write it out!!
  int numF= pMeth_nonCG["F"]["num"];
  int numR= pMeth_nonCG["R"]["num"];
  //std::printf( "%i %i\n", numF, numR);
  // if( (numF == 0) && (numR == 0)) {
  //   if(CpGstatus){std::remove(CpGfile);}
  //   if(CHHstatus){std::remove(CHHfile);}
  //   if(CHGstatus){std::remove(CHGfile);}
  //   // #################
  //   //die("\nnot enough alignments that pass coverage and phred score thresholds to calculate conversion rates\n EXITING....\n\n");
  //   // ################
  // }


  double AvFconvRate  = 0 , AvRconvRate  = 0;
  //int medFconvRate = 0,  medRconvRate = 0;
  if(numF > 0) { AvFconvRate = pMeth_nonCG["F"]["total"]/numF; }
  if(numR > 0) { AvRconvRate = pMeth_nonCG["R"]["total"]/numR; }
  double AvconvRate =(pMeth_nonCG["F"]["total"] + pMeth_nonCG["R"]["total"])/(numF+numR);
  
  //std::vector<double> allesSchon; allesSchon.push_back(pMeth_nonCG["F")); allesSchon.push_back(pMeth_nonCG["R"));
  //if( numF>0 ) { medFconvRate=median(pMeth_nonCG["F"]); }
  //if( numR>0 ) { medRconvRate=median(pMeth_nonCG["R"]); }
  //double medconvRate = median(\@allesSchon);
  
  //int totCs=allesSchon.size();
  int totCs = numF + numR;
  
   
  std::cout 
    <<  "total otherC considered (>95%% C+T): "           << totCs       << "\n"
    <<  "average conversion rate = "                      << AvconvRate  << "\n"
      //"median conversion rate = %.2f\n\n";              //$medconvRate
    
    <<  "total otherC considered (Forward) (>95%% C+T): " << numF        << "\n"
    <<  "average conversion rate (Forward) = "            << AvFconvRate << "\n"
      //"median conversion rate (Forward) = %.2f\n\n"        //$medFconvRate 
    
    <<  "total otherC considered (Reverse) (>95%% C+T): " << numR        << "\n"
    <<  "average conversion rate (Reverse) = "            << AvRconvRate << "\n";
      //"median conversion rate (Reverse) = %.2f\n";       //$medRconvRate
  
  //open (my $hd,">".$prefix."_conversionRate.txt");
  //print $hd $res;
  //close $hd;
  
  return 0;
}



void process_paired_bismark () //std::istream& fh, char* CpGfile, char* CHHfile, char* CHGfile, int &offset, int &mincov, int &minqual, int &nolap, int &paired);
{

printf("Feature is not ready yet.\n");

}

int main() {


  int offset = 33, mincov = 10 , minqual= 20, no_lap = 0, paired = 0;
  
  char buffer[] = "test.fastq_bismark.sorted.min.sam.CpGfile";
  char* empty, *cpgfile =buffer;
  std::ifstream ifs ("test.fastq_bismark.sorted.min.sam" , std::ifstream::in);
  process_sam( std::cin,
   cpgfile,
    empty,
     empty,
      offset,
       mincov,
        minqual,
         no_lap,
          paired);

  return 0;

}



