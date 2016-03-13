#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <deque>
#include <string>
#include <map>
#include <getopt.h>
#include <regex.h>
#include "sam.h"
#include "bgzf.h"
#include "zlib.h"
#include <Rcpp.h>

using namespace Rcpp;


 // ############  Helper functions ################

//  // taken from http://www.kumobius.com/2013/08/c-string-to-int/
// bool String2Int(const std::string& str, int& result)
// {
//     try
//     {
//         std::size_t lastChar;
//         result = std::stoi(str, &lastChar, 10);
//         return lastChar == str.size();
//     }
//     catch (std::invalid_argument&)
//     {
//         return false;
//     }
//     catch (std::out_of_range&)
//     {
//         return false;
//     }
// }


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




// void print_usage() {
//   printf(
//   "usage: ./methcall [options] input_file >outFile.bed"
//   "\n "         
//   "options: \n"
//   " --help     : Print this help message \n"
//   " --read1    : Must be provided at all cases, if given '-' the STDIN will be the input \n"
//   " --type     : one of the following: 'single_sam','paired_sam','single_bismark','paired_bismark' \n"
//   " --nolap    : if given and if the input is paired the overlapping paired reads will be ignored \n"
//   " --minqual  : minquality   (default:20) \n"
//   " --mincov   : min coverage (default:10) \n"
//   " --phred64  : quality scores phred64 scale used otherwise phred33 is the default \n"
//   " --CpG      : output filename for CpG methylation scores (if not specified no file is written out) \n"
//   " --CHH      : output filename for CHH methylation scores (if not specified no file is written out) \n"
//   " --CHG      : output filename for CHG methylation scores (if not specified no file is written out) \n"
//   "           \n"
//   "IMPORTANT: \n"
//   "  Files must be sorted based on chr and start of reads. \n"
//   "  In case of paired-end sam file from bismark, the file still must be sorted in the same way.\n "
//   "\n "
//   );
//   }


// // function to parse command line  arguments
// int get_args(char *&read1, char *&type, int &nolap, int &minqual,
//              int &mincov, int &phred64, char *&CpGfile, char *&CHHfile, char *&CHGfile,
//              int argc, char** argv) {
// 
//   static const struct option longOpts[] = {
//     { "read1",    required_argument,  NULL,     'r' },
//     { "type",     required_argument,  NULL,     't' },
//     { "nolap",    no_argument,        &nolap,    1  },
//     { "minqual",  required_argument,  NULL,     'q' },
//     { "mincov",   required_argument,  NULL,     'c' },
//     { "phred64",  no_argument,        &phred64,  1  },
//     { "CpG",      required_argument,  NULL,     'P' },
//     { "CHH",      required_argument,  NULL,     'H' },
//     { "CHG",      required_argument,  NULL,     'G' },
//     { "help",     no_argument,        NULL,     'h' },
//     { NULL,       no_argument,        NULL,      0  }
//   };
//   
//   static const char *optString = "-r:t:q:c:P:H:G:h";
//   
//   int longIndex = 0;
//   int opt;
// 
//   opt = getopt_long_only(argc,argv, optString, longOpts, &longIndex );
//   
//   while( opt != -1 ) {
//     switch( opt ) {
//       case 'r':   // read in the location
//                   read1 = optarg; 
//                   //std::cout << "set read1 to: " << read1 << std::endl;
//                   break;
//       case 't':   // read in type 
//                   type = optarg;
//                   //std::cout << "set type to: " << type << std::endl;
//                   break;
//       case 'q':   // read minqual and convert to int
//                   minqual = atoi(optarg);
//                   //std::cout << "set minqual to: " << minqual << std::endl;
//                   break;
//       case 'c':   //read mincov and convert to int 
//                   mincov = atoi(optarg);
//                   // std::cout << "set mincov to: " << mincov << std::endl;
//                   break;
//       case 'P':   CpGfile = optarg;
//                   // std::cout << "set CpGfile to: " << CpGfile << std::endl;
//                   break; 
//       case 'H':   CHHfile = optarg;
//                   // std::cout << "set CHHfile to: " << CHHfile << std::endl;
//                   break; 
//       case 'G':   CHGfile = optarg;
//                   // std::cout << "set CHGfile to: " << CHGfile << std::endl;
//                   break;
//       case 'h':   print_usage(argv[0]);
//                   return -1;
//                   break;
//       case ':':   /* missing option argument */
//                   fprintf(stderr, "%s: option `-%c' requires an argument\n",
//                           argv[0], optopt);
//                   return -1;
//                   break;
//       case '?':   /* invalid option */
//                   if ( char(optopt) != 'h') 
//                   { 
//                     fprintf(stderr, "\n%s: option `-%c' is invalid: ignored\n",
//                             argv[0], optopt);
//                   }
//                   break;
//       case 0:
//                 break;     /* getopt_long() set a variable, just keep going */
//       default:    /* You won't actually get here. */
//                   break;
//     }
//       opt = getopt_long_only( argc, argv, optString, longOpts, &longIndex );
//   }
//   
//   return 0;
// }

// check whether required arguments where given and valid
int check_args (const char *read1, const char *type, std::istream *&input, std::ifstream &file) {
     /**
  ###########################################
  # check if required arguments where given #
  ###########################################
  **/

  // check read1 argument
  if( read1==NULL ) 
  { 
      //print_usage();
      Rcpp::stop (" --read1 argument not supplied\n");
      return -1;
  } else {
      //std::cout << read1 << std::endl;
    //string line;
    //reading input from read1
      if( strcmp(read1,"-") ==0 ) {  
        input = &std::cin;
        // tmp = std::cin;
        //while (getline(*tmp,line)) { cout << line << endl; }
      } else {
          file.open(read1);
          if(!file.good()) {Rcpp::stop(" the value of --read1 argument does not point to an existing file\n");return -1;}
          else { 
             input = &file;
            //while (getline(*tmp,line)) { cout << line << endl; } 
          }
        }
    } 

  // check types argument
  std::vector<std::string> types;
  types.push_back("single_sam");
  types.push_back("paired_sam");
  types.push_back("single_bismark");
  types.push_back("paired_bismark");
  types.push_back("bam");

  if( type==NULL) 
  { 
    //print_usage(); 
    Rcpp::stop(" --type argument not supplied\n");
    return -1;
  } else {
      //std::cout << type << std::endl;
    // find returns end of range if element is not found
      if( ( find(types.begin(), types.end(), type)) == types.end()) {
      Rcpp::stop(" --type argument must be one of the following: 'single_sam','paired_sam','single_bismark','paired_bismark','bam' \n");
    }
  }
return 0;
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


// int processnonCGmethHash ( std::map<std::string, std::vector<int> > &nonCGmethHash,std::map<std::string,std::map<std::string, double> >  &CTconvArray, int &mincov) { 
// 
//   // iterate over map/hash
//   for ( std::map<std::string,std::vector<int> >::iterator iter=nonCGmethHash.begin(); iter!=nonCGmethHash.end(); ++iter ) {
//     // get keys and values pair
//     std::string key = iter->first; 
//     std::vector<int> value =  iter->second;
//     // split key to get strand, chr, and locus
//     std::vector<std::string> temp = split(key,'|');
//     std::string strand = temp[0] , chr =  temp[1], loc = temp[2];
//     //int loc = std::atoi(temp[2]); //nextBase; 
//     
//     int noCs= value[0];
//     int noTs= value[1];
//     int noOs= value[2];
//     // double Cperc= 100*noCs/(double)(noTs + noCs + noOs) ;
//     // double Tperc= 100*noTs/(double)(noTs + noCs + noOs) ;
//     if((( noTs + noCs)/(double)(noTs + noCs + noOs) > 0.9) && ((noTs+noCs+noOs)>=mincov) ){
//       //// print to file : "chr.loc   chr   loc   strand   totalbases   %C   %T"
//       //fprintf(out, "%s.%s\t%s\t%s\t%s\t%i\t%.2f\t%.2f\n", 
//               //chr.c_str(), loc.c_str(), chr.c_str(), loc.c_str(), strand.c_str(), noCs+noTs+noOs, Cperc,Tperc); 
//     //}
//     //   push @{$CTconvArray->{$strand}},(($noTs*100)/($noTs+$noCs+$noOs));
//       //CTconvArray[strand]["median"] += ((noTs*100)/(noTs+noCs+noOs));
//     }
//     }
//     return 1;
//     }

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
    CTconvArray[strand]["total"] += ((noTs*100)/(double)(noTs+noCs+noOs));
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


void split_cigar(std::string& source, std::string& regexString, std::vector<std::pair<int,std::string> > &result) {

  size_t maxMatches =  source.length();
  size_t maxGroups = 3;
  regex_t regexCompiled;
  regmatch_t groupArray[maxGroups];
  unsigned int m;
  unsigned int op_len;
  std::string cursor, op;
  
  if (regcomp(&regexCompiled, regexString.c_str(), REG_EXTENDED))
  {
    Rcpp::stop("Could not compile regular expression.\n");
  };
  
  m = 0;
  op_len = 0;
  cursor = source;

  for (m = 0; m < maxMatches; m ++)
  {
    if (regexec(&regexCompiled, cursor.c_str(), maxGroups, groupArray, 0)) break;  // No more matches
    
    // unsigned int g = 0;
    size_t offset = 0;
    
    offset = groupArray[0].rm_eo;
    
    op_len = atoi(cursor.substr(groupArray[1].rm_so,groupArray[1].rm_eo-groupArray[1].rm_so).c_str());
    op = cursor.substr(groupArray[2].rm_so,groupArray[2].rm_eo-groupArray[2].rm_so);
    
//     std::cout << m  << " " <<  cursor 
//               <<  " " << cursor.substr(groupArray[1].rm_so,groupArray[1].rm_eo-groupArray[1].rm_so) 
//               <<  " " << cursor.substr(groupArray[2].rm_so,groupArray[2].rm_eo-groupArray[2].rm_so)
//               << "\n";
    
    result.push_back(std::make_pair(op_len,op));
    cursor.erase(cursor.begin(),cursor.begin()+offset);
  }
  regfree(&regexCompiled);
}

// processes the cigar string and remove and delete elements from mcalls and quality scores
void processCigar ( std::string cigar, std::string &methc, std::string &qual) {
  
  int position = 0, len;
  std::string insertion;
  std::string ops ("MIDS"); // allowed operations
  
  std::vector<std::pair<int,std::string> > cigar_split; // Cigar string is splitted into its single operations
  std::string regexString = "^([0-9]+)([MIDSNHP])";       // -> each pair consists of number and type of operation
  split_cigar(cigar,regexString, cigar_split);          // can be accessed via .first (number) and .second (op)

  std::deque<int> insPos; // location of the insertions
  std::deque<int> insLen; // location of the insert lengths
  
  std::pair<int,std::string> cigar_part;
  
  while (!cigar_split.empty()){
    // if operation is allowed 
    if( ops.find( cigar_split.front().second)!=std::string::npos ) {
      cigar_part =cigar_split.front();
      if (cigar_part.second =="M") { // count matches
        position += cigar_part.first;
      } 
      else if (cigar_part.second == "I") { // count inserts
        len = cigar_part.first;
        insertion = std::string ( len ,'-');
        insPos.push_front(position); 
        insLen.push_front(len); 
        
        position += len;
      } 
      else if (cigar_part.second=="D") { // count deletions
        len = cigar_part.first;
        insertion = std::string(len, '.');
        methc.insert(position, insertion);
        qual.insert(position, insertion);
        
        position += len;
      } 
      else if (cigar_part.second=="S") { 
        //#############################
        Rcpp::stop( "Not ready for this!\n");   
        // ###########################

      }
      // erase the current element
      cigar_split.erase(cigar_split.begin());
    } else {
        // #############################
        Rcpp::stop("Unexpected cigar: "+cigar_part.second+"\n");   
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

// processed sam file without header
int process_sam ( std::istream *fh, std::string &CpGfile, std::string &CHHfile, std::string &CHGfile, int &offset, int &mincov, int &minqual, int nolap, int paired) {

  
// check the file status produce flags 
  int CpGstatus = 0, CHHstatus = 0, CHGstatus = 0;
  FILE *out, *CHHout, *CHGout;
  
  
  if( !(CpGfile.empty()  )) {
    out = fopen(CpGfile.c_str(),"w");
    fprintf (out, "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n");
    CpGstatus=1;
  }
  if( !(CHHfile.empty())) {
    CHHout = fopen(CHHfile.c_str(),"w");
    fprintf (CHHout, "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n");
    CHHstatus=1;
  }
  if( !(CHGfile.empty())) {
    CHGout = fopen(CHGfile.c_str(),"w");
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
  
  while(std::getline(*fh, line))
  {
    
    //std::cout << line << std::endl;
    if(line.find("Bismark") != std::string::npos )  {std::getline(*fh, line);}  // step over the header line
    if(line[0]=='@')                                {std::getline(*fh, line);} // step over the header line
    /** example paired-end reads in SAM format (2 consecutive lines)
    # 1_R1/1	67	5	103172224	255	40M	=	103172417	233	AATATTTTTTTTATTTTAAAATGTGTATTGATTTAAATTT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	NM:i:4	XX:Z:4T1T24TT7	XM:Z:....h.h........................hh.......	XR:Z:CT	XG:Z:CT
    # 1_R1/2	131	5	103172417	255	40M	=	103172224	-233	TATTTTTTTTTAGAGTATTTTTTAATGGTTATTAGATTTT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	NM:i:6	XX:Z:T5T1T9T9T7T3	XM:Z:h.....h.h.........h.........h.......h...	XR:Z:GA	XG:Z:CT
    # HWI-ST986_0098:1:1101:18264:11272#0/1	0	chr1	497	255	50M	*	0	0	TGGGTTTGATTTGAGGAGAATTGTGTTTCGTTTTTAGAGTATTATCGAAA	CCCFFFFFHHHHHJJJIIJJJJJHHIIJIJHIJJJJJGIDHIJJJIIHJI	NM:i:13	XX:Z:C4C3CC9C4C1C2CC2C6CC1C5	XM:Z:z....x...hx.........x....h.xZ.hh..x......hh.xZ....	XR:Z:CT	XG:Z:CT
    **/




    std::vector<std::string> cols = split(line,'\t');
    int start                     = atoi(cols[3].c_str()); 
    // if(!String2Int(cols[3],start)) { Rcpp::stop("Error from String2Int") ;return -1;}
    int end                       = start + cols[9].length()-1;
    std::string chr               = cols[2];
    std::string cigar             = cols[5];
    std::string methc             = cols[13]; 
    methc.erase(methc.begin(),methc.begin()+ 5 ); //  remove "XM:Z:"
    std::string qual              = cols[10];
    std::string mrnm              = cols[6];
    int mpos                      = atoi(cols[7].c_str());
    // int isize                     = atoi(cols[8].c_str());
//     if(!String2Int(cols[7],mpos)) { Rcpp::stop("Error from String2Int");return -1;}
//     if(!String2Int(cols[8],isize)) { Rcpp::stop("Error from String2Int");return -1;}
    
    
    // process cigar string to get indels
    // if search finds nothing string::npos is returned
    if( (cigar.find("D")!=std::string::npos) || (cigar.find("I")!=std::string::npos)) {
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
        Rcpp::stop(  "The sam file is not sorted properly; "
                     "you can sort the file in unix-like machines using:\n" 
                     " grep -v \\'^[[:space:]]*\\@\\' test.sam | sort -k3,3 -k4,4n  > test.sorted.sam \n");
        return -1;
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
  Rcpp::Rcout <<  numF << " " <<  numR << std::endl;
  if( (numF == 0) && (numR == 0)) {
    if(CpGstatus){std::remove(CpGfile.c_str());}
    if(CHHstatus){std::remove(CHHfile.c_str());}
    if(CHGstatus){std::remove(CHGfile.c_str());}
    // #################
    Rcpp::stop("\nnot enough alignments that pass coverage and phred score thresholds to calculate conversion rates\n EXITING....\n\n");
    // ################
    return -1;
  }


  double AvFconvRate  = 0 , AvRconvRate  = 0;
  //int medFconvRate = 0,  medRconvRate = 0;
  if(numF > 0) { AvFconvRate = pMeth_nonCG["F"]["total"]/(double)numF; }
  if(numR > 0) { AvRconvRate = pMeth_nonCG["R"]["total"]/(double)numR; }
  double AvconvRate =(pMeth_nonCG["F"]["total"] + pMeth_nonCG["R"]["total"])/(numF+numR);
  
  //std::vector<double> allesSchon; allesSchon.push_back(pMeth_nonCG["F")); allesSchon.push_back(pMeth_nonCG["R"));
  //if( numF>0 ) { medFconvRate=median(pMeth_nonCG["F"]); }
  //if( numR>0 ) { medRconvRate=median(pMeth_nonCG["R"]); }
  //double medconvRate = median(\@allesSchon);
  
  //int totCs=allesSchon.size();
  int totCs = numF + numR;

   
  Rcpp::Rcout 
    <<  "total otherC considered (>95%% C+T): "           <<   std::setprecision(14) << totCs       << "\n"
    <<  "average conversion rate = "                      <<   std::setprecision(14)<< AvconvRate  << "\n"
      //"median conversion rate = %.2f\n\n";              //$medconvRate
    
    <<  "total otherC considered (Forward) (>95%% C+T): " <<   std::setprecision(14)<< numF        << "\n"
    <<  "average conversion rate (Forward) = "            <<   std::setprecision(14)<< AvFconvRate << "\n"
      //"median conversion rate (Forward) = %.2f\n\n"        //$medFconvRate 
    
    <<  "total otherC considered (Reverse) (>95%% C+T): " <<   std::setprecision(14)<< numR        << "\n"
    <<  "average conversion rate (Reverse) = "            <<   std::setprecision(14)<< AvRconvRate << "\n";
  
  //open (my $hd,">".$prefix."_conversionRate.txt");
  //print $hd $res;
  //close $hd;
  
  return 0;
  
}


// processed  bam/sam file !! with header !!
int process_bam ( std::string &input, std::string &CpGfile, std::string &CHHfile, std::string &CHGfile, int &offset, int &mincov, int &minqual, int nolap) {

  // intialize hts objects which can refer to bam or sam
  htsFile *in = NULL;
  bam1_t *b = NULL;
  bam_hdr_t *header = NULL;

  in = sam_open(input.c_str(),"r");
  if (in==NULL) {Rcpp::stop("fail to open sam/bam file: " + input );}
   
  header = sam_hdr_read(in);
  if ( (header ==NULL) || (header->l_text == 0)) { 
    return 2;
    //Rcpp::stop("fail to open bam header, is there any?\n");
    }

  
// check the file status produce flags 
  int CpGstatus = 0, CHHstatus = 0, CHGstatus = 0;
  FILE *out, *CHHout, *CHGout;
  

  if( !(CpGfile.empty()  )) {
    out = fopen(CpGfile.c_str(),"w");
    fprintf (out, "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n");
    CpGstatus=1;
  }
  if( !(CHHfile.empty())) {
    CHHout = fopen(CHHfile.c_str(),"w");
    fprintf (CHHout, "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n");
    CHHstatus=1;
  }
  if( !(CHGfile.empty())) {
    CHGout = fopen(CHGfile.c_str(),"w");
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

  // initialize bam  
  b = bam_init1();
  //stop("ende");
  // again sam_read1 reads bam or sam files
  while ( sam_read1(in,header,b) >=0 ) 
  {

    
    /** example paired-end reads in SAM format (2 consecutive lines)
    # 1_R1/1  67  5 103172224 255 40M = 103172417 233 AATATTTTTTTTATTTTAAAATGTGTATTGATTTAAATTT  IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII  NM:i:4  XX:Z:4T1T24TT7  XM:Z:....h.h........................hh....... XR:Z:CT XG:Z:CT
    # 1_R1/2  131 5 103172417 255 40M = 103172224 -233  TATTTTTTTTTAGAGTATTTTTTAATGGTTATTAGATTTT  IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII  NM:i:6  XX:Z:T5T1T9T9T7T3 XM:Z:h.....h.h.........h.........h.......h... XR:Z:GA XG:Z:CT
    # HWI-ST986_0098:1:1101:18264:11272#0/1 0 chr1  497 255 50M * 0 0 TGGGTTTGATTTGAGGAGAATTGTGTTTCGTTTTTAGAGTATTATCGAAA  CCCFFFFFHHHHHJJJIIJJJJJHHIIJIJHIJJJJJGIDHIJJJIIHJI  NM:i:13 XX:Z:C4C3CC9C4C1C2CC2C6CC1C5  XM:Z:z....x...hx.........x....h.xZ.hh..x......hh.xZ.... XR:Z:CT XG:Z:CT
    **/


    //get data from bam1_t struct 
    // this link provides some explanations in the definitions
    // https://github.com/samtools/htslib/blob/develop/htslib/sam.h

  int32_t pos = b->core.pos,        // mapping pos aka start 
          len = b->core.l_qseq,     // length of query
          chrom = b->core.tid,      // chromosome id defined by bam_hdr_t, might differ from original names. Compare with ordering header in header!
          mtid = b ->core.mtid,     // chromosome id of next read in template
          mposi = b->core.mpos;     // 0-based leftmost coordinate of next read in template
  
  uint32_t *cigar_pointer = bam_get_cigar(b), // pointer to cigar array
          len_cigar = b->core.n_cigar;        // number of cigar operations

  // read methylation call string
  char *meth = (char*) bam_aux_get(b,"XM");   
  
  // initialize buffers for sequence, qual and cigar string
  std::string seq, qual;
  char cigar_buffer [len];
  
  int i, n;


  // parse the first cigar operation
  n = std::sprintf(cigar_buffer,"%i%c", bam_cigar_oplen(cigar_pointer[0]), bam_cigar_opchr(cigar_pointer[0])) ;
  // if further operations remain, append them to cigar_buffer
  if (len_cigar >= 2 ) { 
    for (i = 1; i < len_cigar; i++  ) {
      char buffer [n];
      n = std::sprintf(buffer,"%i%c", bam_cigar_oplen(cigar_pointer[i]), bam_cigar_opchr(cigar_pointer[i])) ;
      strcat(cigar_buffer,buffer);
    }
  }

  for (i = 0; i < len; i++  ) { qual += bam_get_qual(b)[i]+offset;}


    int start = pos;                 
    int end                       = start + len;
    std::string chr               = header->target_name[chrom]; //get the "real" chromosome id
    std::string cigar             = cigar_buffer ;
    std::string methc             = meth; 
    methc.erase(methc.begin()); //  remove leading "Z" from "XM:Z:"
    // std::string qual              = cols[10];
    std::string mrnm              = header->target_name[mtid];  //get the "real" mate id
    int mpos                      = mposi;
    int paired = (int) ((b)->core.flag&BAM_FPAIRED) ;
    
    // process cigar string to get indels
    // if search finds nothing string::npos is returned
    if( (cigar.find("D")!=std::string::npos) || (cigar.find("I")!=std::string::npos)) {
      processCigar( cigar, methc, qual);
    }
    std::string mcalls = methc; // get the bismark methylation calls
    std::string quals  = qual;  // get the quality scores
    int slen   = mcalls.length(); // aligned sequence length
    
    // get strand
    char strand;
    std::string xr_tag = (char*) bam_aux_get(b, "XR"); //get "XR:Z:CT" from bam/sam, 
    std::string xg_tag = (char*) bam_aux_get(b, "XG");   // but returns "ZCT" as string

    if( (xr_tag == "ZCT") && (xg_tag == "ZCT") ) {strand='+';} // original top strand
    else if( (xr_tag == "ZCT") && (xg_tag == "ZGA") ) {strand='-';} // original bottom strand
    else if( (xr_tag == "ZGA") && (xg_tag == "ZCT") ) {strand='+';} // complementary to original top strand, bismark says - strand to this
    else if( (xr_tag == "ZGA") && (xg_tag == "ZGA") ) {strand='-';} // complementary to original bottom strand, bismark says + strand to this
    
    // if there is no_overlap trim the mcalls and $quals
    // adjust the start
    if( nolap && ( ( mrnm == chr) && paired ) ){
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
        Rcpp::stop ( "The sam file is not sorted properly; "
                     "you can sort the file in unix-like machines using:\n" 
                     " grep -v \\'^[[:space:]]*\\@\\' test.sam | sort -k3,3 -k4,4n  > test.sorted.sam \n");
        return -1;
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
  Rcpp::Rcout <<  numF << " " <<  numR << std::endl;
  if( (numF == 0) && (numR == 0)) {
    if(CpGstatus){std::remove(CpGfile.c_str());}
    if(CHHstatus){std::remove(CHHfile.c_str());}
    if(CHGstatus){std::remove(CHGfile.c_str());}
    // #################
    Rcpp::stop("\nnot enough alignments that pass coverage and phred score thresholds to calculate conversion rates\n EXITING....\n\n");
    // ################
    return -1;
  }


  double AvFconvRate  = 0 , AvRconvRate  = 0;
  //int medFconvRate = 0,  medRconvRate = 0;
  if(numF > 0) { AvFconvRate = pMeth_nonCG["F"]["total"]/(double)numF; }
  if(numR > 0) { AvRconvRate = pMeth_nonCG["R"]["total"]/(double)numR; }
  double AvconvRate =(pMeth_nonCG["F"]["total"] + pMeth_nonCG["R"]["total"])/(numF+numR);
  
  //std::vector<double> allesSchon; allesSchon.push_back(pMeth_nonCG["F")); allesSchon.push_back(pMeth_nonCG["R"));
  //if( numF>0 ) { medFconvRate=median(pMeth_nonCG["F"]); }
  //if( numR>0 ) { medRconvRate=median(pMeth_nonCG["R"]); }
  //double medconvRate = median(\@allesSchon);
  
  //int totCs=allesSchon.size();
  int totCs = numF + numR;

   
  Rcpp::Rcout 
    <<  "total otherC considered (>95%% C+T): "           <<   std::setprecision(14) << totCs       << "\n"
    <<  "average conversion rate = "                      <<   std::setprecision(14)<< AvconvRate  << "\n"
      //"median conversion rate = %.2f\n\n";              //$medconvRate
    
    <<  "total otherC considered (Forward) (>95%% C+T): " <<   std::setprecision(14)<< numF        << "\n"
    <<  "average conversion rate (Forward) = "            <<   std::setprecision(14)<< AvFconvRate << "\n"
      //"median conversion rate (Forward) = %.2f\n\n"        //$medFconvRate 
    
    <<  "total otherC considered (Reverse) (>95%% C+T): " <<   std::setprecision(14)<< numR        << "\n"
    <<  "average conversion rate (Reverse) = "            <<   std::setprecision(14)<< AvRconvRate << "\n";
  
  //open (my $hd,">".$prefix."_conversionRate.txt");
  //print $hd $res;
  //close $hd;
  
  bam_destroy1(b);
  bam_hdr_destroy(header);
  sam_close(in);
  return 0;
  
}




// processed the Bismark 'vanilla' file
int process_single_bismark (std::istream *fh, std::string &CpGfile, std::string &CHHfile, std::string &CHGfile, int &offset, int &mincov, int &minqual ) {

  
// check the file status produce flags 
  int CpGstatus = 0, CHHstatus = 0, CHGstatus = 0;
  FILE *out, *CHHout, *CHGout;
  
  
  if( !(CpGfile.empty()  )) {
    out = fopen(CpGfile.c_str(),"w");
    fprintf (out, "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n");
    CpGstatus=1;
  }
  if( !(CHHfile.empty())) {
    CHHout = fopen(CHHfile.c_str(),"w");
    fprintf (CHHout, "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n");
    CHHstatus=1;
  }
  if( !(CHGfile.empty())) {
    CHGout = fopen(CHGfile.c_str(),"w");
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
  
  while(std::getline(*fh, line))
  {
    
    //std::cout << line << std::endl;
    if(line.find("Bismark") != std::string::npos )  {std::getline(*fh, line);}  // step over the header line
    if(line[0]=='@')                                {std::getline(*fh, line);} // step over the header line
    
    std::vector<std::string> cols = split(line,'\t');
    int start                     = atoi(cols[3].c_str());
    int end                       = atoi(cols[4].c_str());                    
//     if(!String2Int(cols[3],start)) { Rcpp::stop("Error from String2Int");return -1;}
//     if(!String2Int(cols[4],end)) { Rcpp::stop( "Error from String2Int");return -1;}
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
        Rcpp::stop( "The sam file is not sorted properly; you can sort the file in unix-like machines using:\n"
                    "grep -v \\'^[[:space:]]*\\@\\' test.sam | sort -k3,3 -k4,4n  > test.sorted.sam \n");
        return -1;
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
  //fclose(out);


  if(CpGstatus) { processCGmethHash(CGmethHash,out,mincov);     std::fclose(out); }
  if(CHHstatus) { processCHmethHash(CHHmethHash,CHHout,mincov); std::fclose(CHHout); }
  if(CHGstatus) { processCHmethHash(CHGmethHash,CHGout,mincov); std::fclose(CHGout); }
  //processnonCGmethHash(nonCGmethHash,pMeth_nonCG,mincov);
  processnonCGmethHash2(nonCGmethHash,pMeth_nonCG,mincov);


  // get the conversion rate and write it out!!
  int numF= pMeth_nonCG["F"]["num"];
  int numR= pMeth_nonCG["R"]["num"];
  Rcpp::Rcout <<  numF << " " <<  numR << std::endl;
  if( (numF == 0) && (numR == 0)) {
    if(CpGstatus){std::remove(CpGfile.c_str());}
    if(CHHstatus){std::remove(CHHfile.c_str());}
    if(CHGstatus){std::remove(CHGfile.c_str());}
    // #################
    Rcpp::stop("\nnot enough alignments that pass coverage and phred score thresholds "
               "to calculate conversion rates\n EXITING....\n\n");
    // ################
    return -1;
  }


  double AvFconvRate  = 0 , AvRconvRate  = 0;
  //int medFconvRate = 0,  medRconvRate = 0;
  if(numF > 0) { AvFconvRate = pMeth_nonCG["F"]["total"]/(double) numF; }
  if(numR > 0) { AvRconvRate = pMeth_nonCG["R"]["total"]/(double) numR; }
  double AvconvRate =(pMeth_nonCG["F"]["total"] + pMeth_nonCG["R"]["total"])/(numF+numR);
  
  //std::vector<double> allesSchon; allesSchon.push_back(pMeth_nonCG["F")); allesSchon.push_back(pMeth_nonCG["R"));
  //if( numF>0 ) { medFconvRate=median(pMeth_nonCG["F"]); }
  //if( numR>0 ) { medRconvRate=median(pMeth_nonCG["R"]); }
  //double medconvRate = median(\@allesSchon);
  
  //int totCs=allesSchon.size();
  int totCs = numF + numR;
  
   
  Rcpp::Rcout 
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

Rcpp::stop("Feature is not ready yet.\n");

}

// int main(int argc,const char **argv){
// 
//   
//   // initialize variables
//   char *read1 = NULL,   *type = NULL, 
//     *CpGfile = NULL, *CHGfile = NULL, 
//     *CHHfile = NULL;
//     int minqual = 20, mincov  = 10;
//     int phred64 = 0,  nolap = 0;
//     int offset = 33;
// 
//   char
//     
//   
//   if( get_args(read1, type, nolap, minqual, mincov, phred64, CpGfile, CHHfile, CHGfile, argc, (char*) argv) != 0 ) {
//     return -1;
//   }

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
void methCall(std::string read1, std::string type="bam", bool nolap=false, int minqual=20,
                int mincov = 10 , bool phred64 = false , std::string CpGfile ="",
                std::string CHHfile ="" , std::string CHGfile = "" ) {

  int offset = 33;
  if (phred64) offset = 64;

  std::istream *input = NULL;
  std::string line;
  std::ifstream file;
 
  check_args(read1.c_str(), type.c_str(), input, file);
  
  int res = 0;
  
  if(!type.empty()) {
    // type is "bam" per default (reads both sam or bam), but if the header is missing, 
    // the file will be treated as paired_sam
    if( type == "bam"){
      res = process_bam(read1, CpGfile, CHHfile, CHGfile, offset, mincov, minqual ,nolap);
      // std::cout << res;
    }
    if( (res==2)  || (type == "paired_sam")){
      process_sam(input, CpGfile,CHHfile,CHGfile,offset,mincov,minqual,nolap,1);
    }
    else if(type == "single_sam"){
      process_sam(input, CpGfile, CHHfile, CHGfile, offset, mincov, minqual ,0,0);
    }
    else if( type == "single_bismark" ){
      process_single_bismark(input, CpGfile,CHHfile,CHGfile,offset,mincov,minqual);
    }
    else if( type =="paired_bismark"){
    Rcpp::stop( "--paired_bismark option NOT IMPLEMENTED! get a paired sam file and used that as input\n") ;
    // return false;
    }

  }
  
  if(file.is_open()) file.close();

}

