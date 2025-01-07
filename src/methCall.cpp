


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
//#include <getopt.h>
#include "htslib/sam.h"
#include "htslib/bgzf.h"
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
std::vector<std::string> &split(const std::string &s, char delim,
                                std::vector<std::string> &elems) {
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


// small function to replace a substring
std::string find_and_replace(std::string source, std::string const& find, std::string const& replace)
{
    for(std::string::size_type i = 0; (i = source.find(find, i)) != std::string::npos;)
    {
        source.replace(i, find.length(), replace);
        i += replace.length();
    }

    return(source);
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
    if((( noTs + noCs)/(double)(noTs + noCs + noOs) > 0.95) && ((noTs+noCs+noOs)>=mincov) ){
    //print join("\t",($chr.".".$loc,$chr,$loc,$strand,$noCs+$noTs+$noOs,$Cperc,$Tperc)  ),"\n"; 
    //push @{$CTconvArray->{$strand}},(($noTs*100)/($noTs+$noCs+$noOs));
    CTconvArray[strand]["num"]++;
    //CTconvArray[strand]["total"] += ((noTs*100)/(noTs+noCs+noOs));
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
// 
//  The methylation call string contains:
//  - a dot . for every position in the BS-read not involving a cytosine,
//  - one of the following letters for the three different cytosine methylation contexts,
//    where (UPPER CASE = METHYLATED, lower case = unmethylated):
// 
//     z - C in CpG context - unmethylated
//     Z - C in CpG context - methylated
//     x - C in CHG context - unmethylated
//     X - C in CHG context - methylated
//     h - C in CHH context - unmethylated
//     H - C in CHH context - methylated
//     u - C in Unknown context (CN or CHN) - unmethylated
//     U - C in Unknown context (CN or CHN) - methylated
//     . - not a C or irrelevant position
//
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
  
  int position = 0, len, i = 0;
  std::string insertion;
  std::string ops ("MIDS"); // allowed operations
  
  
  std::vector<std::pair<int,std::string> > cigar_split; // Cigar string is splitted into its single operations
                                                        // -> each pair consists of number and type of operation
                                                        // can be accessed via .first (number) and .second (op)
  std::string op;
  int oplen, oplen_start=0, op_pos;
  std::string::const_iterator it;
  for (it = cigar.begin(); it < cigar.end() ; ++it) {
    
    if( std::isalpha(*it) )    {
      op = *it;
      op_pos=it - cigar.begin();
      oplen = atoi(cigar.substr(oplen_start,op_pos).c_str());
      oplen_start = op_pos+1;
      
      //std::cout << oplen << " "<<  op << std::endl;
      cigar_split.push_back(std::make_pair(oplen,op));
    }
  }                                                      

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
    for( i=0; ((unsigned) i) < insPos.size(); i++){
      methc.erase(insPos[i],insLen[i]);
      qual.erase(insPos[i],insLen[i]);
    }
  }
  
}



// processed sam file without header
int process_sam ( std::istream *fh, 
                  std::string &CpGfile, 
                  std::string &CHHfile, 
                  std::string &CHGfile, 
                  int &offset, 
                  int &mincov, 
                  int &minqual, 
                  int nolap, 
                  int paired, 
                  size_t &verbosity) {

  
// check the file status produce flags 
  int CpGstatus = 0, CHHstatus = 0, CHGstatus = 0;
  FILE *out, *CHHout, *CHGout;
  FILE *out_stats, *CHHout_stats, *CHGout_stats;
  std::string tempString ;
  
  
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
  int i = 0; 
  std::string lastChrom="null", chrPre;
  std::vector<std::string> pastChrom; 

  std::string line;  
  
  while(std::getline(*fh, line))
  {
    //check wheter user wnats to interrupt runnning code
    Rcpp::checkUserInterrupt();
    
    //std::cout << line << std::endl;
    if(line.find("Bismark") != std::string::npos )  {std::getline(*fh, line);}  // step over the header line
    while(line[0]=='@')                                {std::getline(*fh, line);} // step over the header line
    /** example paired-end reads in SAM format (2 consecutive lines)
    # 1_R1/1	67	5	103172224	255	40M	=	103172417	233	AATATTTTTTTTATTTTAAAATGTGTATTGATTTAAATTT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	NM:i:4	XX:Z:4T1T24TT7	XM:Z:....h.h........................hh.......	XR:Z:CT	XG:Z:CT
    # 1_R1/2	131	5	103172417	255	40M	=	103172224	-233	TATTTTTTTTTAGAGTATTTTTTAATGGTTATTAGATTTT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	NM:i:6	XX:Z:T5T1T9T9T7T3	XM:Z:h.....h.h.........h.........h.......h...	XR:Z:GA	XG:Z:CT
    # HWI-ST986_0098:1:1101:18264:11272#0/1	0	chr1	497	255	50M	*	0	0	TGGGTTTGATTTGAGGAGAATTGTGTTTCGTTTTTAGAGTATTATCGAAA	CCCFFFFFHHHHHJJJIIJJJJJHHIIJIJHIJJJJJGIDHIJJJIIHJI	NM:i:13	XX:Z:C4C3CC9C4C1C2CC2C6CC1C5	XM:Z:z....x...hx.........x....h.xZ.hh..x......hh.xZ....	XR:Z:CT	XG:Z:CT
    **/




    std::vector<std::string> cols = split(line,'\t');
    // for( i=0; ((unsigned) i)<cols.size(); ++i)
    //   Rcpp::Rcout << cols[i] << ' ';
    // Rcpp::Rcout  << std::endl;  
    int start                     = atoi(cols[3].c_str()); 
    int end                       = start + cols[9].length()-1;
    std::string chr               = cols[2];
    std::string cigar             = cols[5];
    std::string methc             = cols[13];
    if(methc.find("XM:Z:")==std::string::npos) {Rcpp::stop("no methylation tag found.");}
    methc.erase(methc.begin(),methc.begin()+ 5 ); //  remove "XM:Z:"
    std::string qual              = cols[10];
    std::string mrnm              = cols[6];
    int mpos                      = atoi(cols[7].c_str());

    
    
    // process cigar string to get indels
    // if search finds nothing string::npos is returned
    if( (cigar.find("D")!=std::string::npos) || (cigar.find("I")!=std::string::npos)) {
      processCigar( cigar, methc, qual);
    }
    std::string mcalls = methc; // get the bismark methylation calls
    std::string quals  = qual;  // get the quality scores
    int slen   = mcalls.length(); // aligned sequence length
    
    // get strand
    char strand = ' ';
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
        Rcpp::stop(  "The sam file is not sorted properly on positions :\n"
                     "\tchr: %s  pos: %i is followed by chr: %s  pos: %i \n"
                     "You can sort the file in unix-like machines using:\n" 
                     " grep -v \\'^[[:space:]]*\\@\\' test.sam | sort -k3,3 -k4,4n  > test.sorted.sam \n"
                     ,chrPre, startPre, chr,start);
        return -1;
        // ########################
      }
      chrPre=chr;
      startPre=start;
    } else {
      if ( std::find(pastChrom.begin(), pastChrom.end(), chr) != pastChrom.end() ) {
        // ####################### 
        Rcpp::stop(  "The sam file is not sorted properly on chromosomes:\n"
                       "\tchr: %s occured before and after %s.\n"
                       "You can sort the file in unix-like machines using:\n" 
                       " grep -v \\'^[[:space:]]*\\@\\' test.sam | sort -k3,3 -k4,4n  > test.sorted.sam \n"
                       ,chr,chrPre);
        return -1;
        // ########################
      }
      startPre=start;
      chrPre=chr;
      pastChrom.push_back(chr);
      // for( i=0; ((unsigned) i)<pastChrom.size(); ++i)
      //   Rcpp::Rcout << pastChrom[i] << ' ';
      // Rcpp::Rcout  << std::endl;  
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
    for( i=0; ((unsigned) i) < quals.length(); i++) 
    {
      if ( ( ( (int) quals[i] - offset ) <  minqual ) || ( mcalls[i] == '.') ){ continue;}
      std::string key; // initialize the hash key
      if( strand == '+') { key = "F|"+ chr+"|"+std::to_string(static_cast<long long>(start+i)); }
      else { key = "R|"+ chr+"|"+std::to_string(static_cast<long long>(start+i)); }

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
  //Rcpp::Rcout <<  numF << " " <<  numR << std::endl;
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

  std::stringstream sout;
  
  sout   <<   "Conversion Statistics:\n\n" 
   <<  std::setprecision(14) 
   <<  "total otherC considered (>95% C+T): "            <<   totCs       << "\n"
   <<  "average conversion rate = "                      <<   AvconvRate  << "\n"
   
   <<  "total otherC considered (Forward) (>95% C+T): "  <<   numF        << "\n"
   <<  "average conversion rate (Forward) = "            <<   AvFconvRate << "\n"
   
   <<  "total otherC considered (Reverse) (>95% C+T): "  <<   numR        << "\n"
   <<  "average conversion rate (Reverse) = "            <<   AvRconvRate << "\n";


  if(verbosity > 0 ) Rcpp::Rcout << sout.str() << std::endl;
  

  if( !(CpGfile.empty()  )) {
    tempString = CpGfile; 
    tempString.append("_conversionStats.txt");
    out_stats = fopen(find_and_replace( tempString, ".txt_", "_").c_str(),"w");
    fprintf (out_stats, "%s", sout.str().c_str());
    std::fclose(out_stats);
  }
  if( !(CHHfile.empty())) {
    tempString = CHHfile; 
    tempString.append("_conversionStats.txt");
    CHHout_stats = fopen(find_and_replace( tempString, ".txt_", "_").c_str(),"w");
    fprintf (CHHout_stats, "%s", sout.str().c_str());
    std::fclose(CHHout_stats);
  }
  if( !(CHGfile.empty())) {
    tempString = CHGfile; 
    tempString.append("_conversionStats.txt");
    CHGout_stats = fopen(find_and_replace( tempString, ".txt_", "_").c_str(),"w");
    fprintf (CHGout_stats, "%s", sout.str().c_str());
    std::fclose(CHGout_stats);
  }


  return 0;
  
}



// processed  single-end bam/sam file !! with header !!
int process_bam ( std::string &input, 
                  std::string &CpGfile, 
                  std::string &CHHfile, 
                  std::string &CHGfile, 
                  int &offset, 
                  int &mincov, 
                  int &minqual, 
                  int nolap, 
                  size_t &verbosity) {

  // intialize hts objects which can refer to bam or sam
  htsFile *in = NULL;
  bam1_t *b = NULL;
  bam_hdr_t *header = NULL;

  if (verbosity > 1)
    Rcpp::Rcout << "\n"
                << "Processing " << input << "\n\n"
                << std::endl;
  in = sam_open(input.c_str(),"r");
  if (in==NULL) {Rcpp::stop("fail to open sam/bam file: " + input );}
   
  header = sam_hdr_read(in);
  if ( (header ==NULL) || (header->l_text == 0)) { 
    if(verbosity > 1 ) Rcpp::Rcout << "Failed to read header, falling back." << std::endl ; 
    return 2;
    }

  
// check the file status produce flags 
  int CpGstatus = 0, CHHstatus = 0, CHGstatus = 0;
  FILE *out, *CHHout, *CHGout;
  FILE *out_stats, *CHHout_stats, *CHGout_stats;
  std::string tempString;
  

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
  
  // initialize maps to store intermediate data
  std::map<std::string, std::vector<int> > CGmethHash; 
  std::map<std::string, std::vector<int> > nonCGmethHash; 
  std::map<std::string, std::vector<int> > CHHmethHash;
  std::map<std::string, std::vector<int> > CHGmethHash;
  std::map<std::string,std::map<std::string, double> > pMeth_nonCG;

  
  int lastPos  =-1, startPre = -1 ;
  std::string lastChrom="null", chrPre;
  std::vector<std::string> pastChrom; 

  // initialize bam  
  b = bam_init1();
  //stop("ende");
  // again sam_read1 reads bam or sam files
  while ( sam_read1(in,header,b) >=0 ) 
  {

    //check wheter user wnats to interrupt runnning code
    Rcpp::checkUserInterrupt();
    
    /** example paired-end reads in SAM format (2 consecutive lines)
    # 1_R1/1  67  5 103172224 255 40M = 103172417 233 AATATTTTTTTTATTTTAAAATGTGTATTGATTTAAATTT  IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII  NM:i:4  XX:Z:4T1T24TT7  XM:Z:....h.h........................hh....... XR:Z:CT XG:Z:CT
    # 1_R1/2  131 5 103172417 255 40M = 103172224 -233  TATTTTTTTTTAGAGTATTTTTTAATGGTTATTAGATTTT  IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII  NM:i:6  XX:Z:T5T1T9T9T7T3 XM:Z:h.....h.h.........h.........h.......h... XR:Z:GA XG:Z:CT
    # HWI-ST986_0098:1:1101:18264:11272#0/1 0 chr1  497 255 50M * 0 0 TGGGTTTGATTTGAGGAGAATTGTGTTTCGTTTTTAGAGTATTATCGAAA  CCCFFFFFHHHHHJJJIIJJJJJHHIIJIJHIJJJJJGIDHIJJJIIHJI  NM:i:13 XX:Z:C4C3CC9C4C1C2CC2C6CC1C5  XM:Z:z....x...hx.........x....h.xZ.hh..x......hh.xZ.... XR:Z:CT XG:Z:CT
    **/


    //get data from bam1_t struct 
    // this link provides some explanations in the definitions
    // https://github.com/samtools/htslib/blob/develop/htslib/sam.h

    char *qname = (char *)bam_get_qname(b);
    if (verbosity > 1)
      std::printf("read name: %s\n", qname);

    int32_t pos = b->core.pos, // 0-based leftmost coordinate
        len = b->core.l_qseq,  // length of query
        chrom = b->core.tid,   // chromosome id defined by bam_hdr_t, might differ from original names. Compare with ordering header in header!
        mtid = b->core.mtid,   // chromosome id of next read in template
        mpos = b->core.mpos;   // 0-based leftmost coordinate of next read in template

    // change pos and mposi to 1-based coordinates to get same results as for sam input
    pos = pos + 1;
    mpos = (int)mpos + 1;

    uint32_t *cigar_pointer = bam_get_cigar(b), // pointer to cigar array
        len_cigar = b->core.n_cigar;            // number of cigar operations

    // read methylation call string
    char *meth = (char *)bam_aux_get(b, "XM");
    if (meth == NULL)
    {
      Rcpp::stop("no methylation tag found for bam file " + input);
    }

    if (verbosity > 1)
      std::cout << "pos: " << pos << " len: " << len << " chrom: " << chrom << " mtid: " << mtid << " mpos: " << mpos << " len_cigar: " << len_cigar << std::endl;

    // initialize buffers for sequence, qual and cigar string
    std::string seq, qual;
    std::string cigar;

    // initialize counter and cigar-buffer-storage
    int i = 0;
    int c;
    std::string cigar_buffer;

    // parse the cigar operations
    for (i = 0; i < (int)len_cigar; i++)
    {
      // format the cigar operations as string and save the string length as c
      c = std::snprintf(cigar_buffer.data(), cigar_buffer.size(), "%i%c", bam_cigar_oplen(cigar_pointer[i]), bam_cigar_opchr(cigar_pointer[i]));
      if (verbosity > 1)
      {
        std::printf("cigar buffer length: %i\n", c);
        std::printf("cigar operation length: %i\n", bam_cigar_oplen(cigar_pointer[i]));
        std::printf("cigar operation performed: %c\n", bam_cigar_opchr(cigar_pointer[i]));
      }

      // place each operation as character into cigar_buffer
      std::snprintf(cigar_buffer.data(), c + 1, "%i%c", bam_cigar_oplen(cigar_pointer[i]), bam_cigar_opchr(cigar_pointer[i]));
      // append to cigar string
      cigar += cigar_buffer.c_str();

      // print current cigar_buffer
      if (verbosity > 1)
      {
        std::printf("cigar operation: %i\n", i);
        std::printf("cigar_buffer: %s\n", cigar_buffer.c_str());
        std::printf("cigar string: %s\n", cigar.c_str());
      }
    }

    for (i = 0; i < len; i++)
    {
      qual += bam_get_qual(b)[i] + offset;
    }

    int start = pos;                 
    int end                       = start + len;
    std::string chr = header->target_name[chrom]; // get the "real" chromosome id
    std::string methc             = meth; 
    methc.erase(methc.begin()); //  remove leading "Z" from "XM:Z:"
    // std::string qual              = cols[10];

    // print out some information
    if (verbosity > 1)
      std::cout << "start: " << start << " end: " << end << " chr: " << chr << " cigar: " << cigar << " meth: " << methc << std::endl;

    // // check BAM flag if read was paired in sequencing
    // int paired = (int) ((b)->core.flag&BAM_FPAIRED) ;
    
    // process cigar string to get indels
    // if search finds nothing string::npos is returned
    if( (cigar.find("D")!=std::string::npos) || (cigar.find("I")!=std::string::npos)) {
      processCigar( cigar, methc, qual);
    }
    std::string mcalls = methc; // get the bismark methylation calls
    std::string quals  = qual;  // get the quality scores
    int slen   = mcalls.length(); // aligned sequence length
    
    // get strand
    char strand = ' ';
    std::string xr_tag = (char*) bam_aux_get(b, "XR"); //get "XR:Z:CT" from bam/sam, 
    std::string xg_tag = (char*) bam_aux_get(b, "XG");   // but returns "ZCT" as string

    if( (xr_tag == "ZCT") && (xg_tag == "ZCT") ) {strand='+';} // original top strand
    else if( (xr_tag == "ZCT") && (xg_tag == "ZGA") ) {strand='-';} // original bottom strand
    else if( (xr_tag == "ZGA") && (xg_tag == "ZCT") ) {strand='+';} // complementary to original top strand, bismark says - strand to this
    else if( (xr_tag == "ZGA") && (xg_tag == "ZGA") ) {strand='-';} // complementary to original bottom strand, bismark says + strand to this
    
    
    // check wether read is proper paired (both mapped)
    int proper_paired = (int) ((b)->core.flag&BAM_FPROPER_PAIR);
    
    // if is proper pair and no_overlap is set
    // -> trim the mcalls and $quals
    // -> adjust the start
    if( nolap && proper_paired  ){
      std::string mrnm = header->target_name[mtid];  //get the mate id
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
        Rcpp::stop(  "The sam file is not sorted properly on positions :\n"
                     "\tchr: %s  pos: %i is followed by chr: %s  pos: %i \n"
                     "You can sort the file in unix-like machines using:\n" 
                     " grep -v \\'^[[:space:]]*\\@\\' test.sam | sort -k3,3 -k4,4n  > test.sorted.sam \n"
                     ,chrPre, startPre, chr,start);
        return -1;
        // ########################
      }
      chrPre=chr;
      startPre=start;
    } else {
      if ( std::find(pastChrom.begin(), pastChrom.end(), chr) != pastChrom.end() ) {
        // ####################### 
        Rcpp::stop(  "The sam file is not sorted properly on chromosomes:\n"
                       "\tchr: %s occured before and after %s.\n"
                       "You can sort the file in unix-like machines using:\n" 
                       " grep -v \\'^[[:space:]]*\\@\\' test.sam | sort -k3,3 -k4,4n  > test.sorted.sam \n"
                       ,chr,chrPre);
        return -1;
        // ########################
      }
      startPre=start;
      chrPre=chr;
      pastChrom.push_back(chr);
      // for( i=0; ((unsigned) i)<pastChrom.size(); ++i)
      //   Rcpp::Rcout << pastChrom[i] << ' ';
      // Rcpp::Rcout  << std::endl;  
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
    for(  i=0; ((unsigned) i) < quals.length(); i++) 
    {
      if ( ( ( (int) quals[i] - offset ) <  minqual ) || ( mcalls[i] == '.') ){ continue;}
      std::string key; // initialize the hash key
      if( strand == '+') { key = "F|"+ chr+"|"+std::to_string(static_cast<long long>(start+i)); }
      else { key = "R|"+ chr+"|"+std::to_string(static_cast<long long>(start+i)); }

      process_call_string(mcalls,i,key, CGmethHash, nonCGmethHash, CHHmethHash, CHGmethHash);
      
      
    }
    lastPos=end;
    lastChrom=chr;

    if (verbosity > 1)
      Rcpp::Rcout << "read done\n"
                  << std::endl;
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
  //Rcpp::Rcout <<  numF << " " <<  numR << std::endl;
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
  
  int totCs = numF + numR;
   
  std::stringstream sout;
  
  sout   <<   "Conversion Statistics:\n\n" 
   <<  std::setprecision(14) 
   <<  "total otherC considered (>95% C+T): "            <<   totCs       << "\n"
   <<  "average conversion rate = "                      <<   AvconvRate  << "\n"
   
   <<  "total otherC considered (Forward) (>95% C+T): "  <<   numF        << "\n"
   <<  "average conversion rate (Forward) = "            <<   AvFconvRate << "\n"
   
   <<  "total otherC considered (Reverse) (>95% C+T): "  <<   numR        << "\n"
   <<  "average conversion rate (Reverse) = "            <<   AvRconvRate << "\n";


  if(verbosity > 0 ) Rcpp::Rcout << sout.str() << std::endl;
  

  if( !(CpGfile.empty()  )) {
    tempString = CpGfile; 
    tempString.append("_conversionStats.txt");
    out_stats = fopen(find_and_replace( tempString, ".txt_", "_").c_str(),"w");
    fprintf (out_stats, "%s", sout.str().c_str());
    std::fclose(out_stats);
  }
  if( !(CHHfile.empty())) {
    tempString = CHHfile; 
    tempString.append("_conversionStats.txt");
    CHHout_stats = fopen(find_and_replace( tempString, ".txt_", "_").c_str(),"w");
    fprintf (CHHout_stats, "%s", sout.str().c_str());
    std::fclose(CHHout_stats);
  }
  if( !(CHGfile.empty())) {
    tempString = CHGfile; 
    tempString.append("_conversionStats.txt");
    CHGout_stats = fopen(find_and_replace( tempString, ".txt_", "_").c_str(),"w");
    fprintf (CHGout_stats, "%s", sout.str().c_str());
    std::fclose(CHGout_stats);
  }

  
  bam_destroy1(b);
  bam_hdr_destroy(header);
  sam_close(in);
  return 0;
  
}




// processed the Bismark 'vanilla' file
int process_single_bismark (std::istream *fh, 
                            std::string &CpGfile, 
                            std::string &CHHfile, 
                            std::string &CHGfile, 
                            int &offset, 
                            int &mincov, 
                            int &minqual, 
                            size_t &verbosity) {

  
// check the file status produce flags 
  int CpGstatus = 0, CHHstatus = 0, CHGstatus = 0;
  FILE *out, *CHHout, *CHGout;
  FILE *out_stats, *CHHout_stats, *CHGout_stats;
  std::string tempString;
  
  
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
  
  
  
  int i = 0;
  int lastPos  =-1, startPre = -1 ;
  std::string lastChrom="null", chrPre;
  std::vector<std::string> pastChrom; 

  std::string line;  
  
  while(std::getline(*fh, line))
  {
    
    //check wheter user wnats to interrupt runnning code
    Rcpp::checkUserInterrupt();
    
    //std::cout << line << std::endl;
    if(line.find("Bismark") != std::string::npos )  {std::getline(*fh, line);}  // step over the header line
    while(line[0]=='@')                                {std::getline(*fh, line);} // step over the header line
    
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
        Rcpp::stop(  "The sam file is not sorted properly on positions :\n"
                     "\tchr: %s  pos: %i is followed by chr: %s  pos: %i \n"
                     "You can sort the file in unix-like machines using:\n" 
                     " grep -v \\'^[[:space:]]*\\@\\' test.sam | sort -k3,3 -k4,4n  > test.sorted.sam \n"
                     ,chrPre, startPre, chr,start);
        return -1;
        // ########################
      }
      chrPre=chr;
      startPre=start;
    } else {
      if ( std::find(pastChrom.begin(), pastChrom.end(), chr) != pastChrom.end() ) {
        // ####################### 
        Rcpp::stop(  "The sam file is not sorted properly on chromosomes:\n"
                       "\tchr: %s occured before and after %s.\n"
                       "You can sort the file in unix-like machines using:\n" 
                       " grep -v \\'^[[:space:]]*\\@\\' test.sam | sort -k3,3 -k4,4n  > test.sorted.sam \n"
                       ,chr,chrPre);
        return -1;
        // ########################
      }
      startPre=start;
      chrPre=chr;
      pastChrom.push_back(chr);
      // for( i=0; ((unsigned) i)<pastChrom.size(); ++i)
      //   Rcpp::Rcout << pastChrom[i] << ' ';
      // Rcpp::Rcout  << std::endl;  
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
    for( i=0; ((unsigned) i) < quals.length(); i++) 
    {
      if ( ( ( (int) quals[i] - offset ) <  minqual ) || ( mcalls[i] == '.') ){ continue;}
      //if last base is a C and it is a part of CCGG motif, don't call for meth
      if( ( (gbases[i] == 'C') && ( ((unsigned) i) == quals.length()) ) && ( gbases.substr(i-1,4) == "CCGG" ) ) { continue;} 
      std::string key; // initilaize the hash key
      if( strand == '+') { key = "F|"+ chr+"|"+std::to_string(static_cast<long long>(start+i)); }
      else { key = "R|"+ chr+"|"+std::to_string(static_cast<long long>(start+i)); }
      
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
  // Rcpp::Rcout <<  numF << " " <<  numR << std::endl;
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
  
   
  std::stringstream sout;
  
  sout   <<   "Conversion Statistics:\n\n" 
   <<  std::setprecision(14) 
   <<  "total otherC considered (>95% C+T): "            <<   totCs       << "\n"
   <<  "average conversion rate = "                      <<   AvconvRate  << "\n"
   
   <<  "total otherC considered (Forward) (>95% C+T): "  <<   numF        << "\n"
   <<  "average conversion rate (Forward) = "            <<   AvFconvRate << "\n"
   
   <<  "total otherC considered (Reverse) (>95% C+T): "  <<   numR        << "\n"
   <<  "average conversion rate (Reverse) = "            <<   AvRconvRate << "\n";


  if(verbosity > 0 ) Rcpp::Rcout << sout.str() << std::endl;
  

  if( !(CpGfile.empty()  )) {
    tempString = CpGfile; 
    tempString.append("_conversionStats.txt");
    out_stats = fopen(find_and_replace( tempString, ".txt_", "_").c_str(),"w");
    fprintf (out_stats, "%s", sout.str().c_str());
    std::fclose(out_stats);
  }
  if( !(CHHfile.empty())) {
    tempString = CHHfile; 
    tempString.append("_conversionStats.txt");
    CHHout_stats = fopen(find_and_replace( tempString, ".txt_", "_").c_str(),"w");
    fprintf (CHHout_stats, "%s", sout.str().c_str());
    std::fclose(CHHout_stats);
  }
  if( !(CHGfile.empty())) {
    tempString = CHGfile; 
    tempString.append("_conversionStats.txt");
    CHGout_stats = fopen(find_and_replace( tempString, ".txt_", "_").c_str(),"w");
    fprintf (CHGout_stats, "%s", sout.str().c_str());
    std::fclose(CHGout_stats);
  }

  
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
                std::string CHHfile ="" , std::string CHGfile = "" , size_t verbosity = 1) {

  int offset = 33;
  if (phred64) offset = 64;

  std::istream *input = NULL;
  //std::string line;
  std::ifstream file;
 
  check_args(read1.c_str(), type.c_str(), input, file);
  
  int res = 0;
  
  if(verbosity > 0 ) Rcpp::Rcout << "Trying to process:\n\t" << read1 << std::endl;
  
  if(!type.empty()) {
    // type is "bam" per default (reads both sam or bam), but if the header is missing, 
    // the file will be treated as paired_sam
    if( type == "bam"){
      if(verbosity > 1 ) Rcpp::Rcout << "Using htslib."  << std::endl;
      res = process_bam(read1, CpGfile, CHHfile, CHGfile, offset, mincov, minqual ,nolap, verbosity);
      // std::cout << res;
    }
    if( (res==2)  || (type == "paired_sam")){
      if(verbosity > 1 ) Rcpp::Rcout << "As paired sam."  << std::endl;
      process_sam(input, CpGfile,CHHfile,CHGfile,offset,mincov,minqual,nolap,1,verbosity);
    }
    else if(type == "single_sam"){
      if(verbosity > 1 ) Rcpp::Rcout << "As single sam."  << std::endl;
      process_sam(input, CpGfile, CHHfile, CHGfile, offset, mincov, minqual ,0,0,verbosity);
    }
    else if( type == "single_bismark" ){
      if(verbosity > 1 ) Rcpp::Rcout << "As single bismark."  << std::endl;
      process_single_bismark(input, CpGfile,CHHfile,CHGfile,offset,mincov,minqual,verbosity);
    }
    else if( type =="paired_bismark"){
    Rcpp::stop( "--paired_bismark option NOT IMPLEMENTED! get a paired sam file and used that as input\n") ;
    // return false;
    }

  }

  if (verbosity > 1)
    Rcpp::Rcout << "Done.\n"
                << std::endl;

  if(file.is_open()) file.close();

}
