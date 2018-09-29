//
//  main.cpp
//  GenDuplexGroup
//
//  Created by Lee on 2017/7/3.
//  Copyright © 2017年 Lee. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <getopt.h>

#include "TimeFunc.hpp"
#include "Files.hpp"
#include "environment.hpp"
#include "DuplexGroup.hpp"

using namespace std;


//void genDuplexGroup(string read_pair_file, string duplex_group_file, int OVERLAP=5, bool multiDG=false);
//void collapseDuplexGroup(string sorted_duplex_group_file, string collapsed_dg_file, int maxGap=10, int maxTotal=30);


#define PARIS_VERSION "1.0"
#define BUILD_HOST "Li Pan"
#define BUILD_TIME "2017-7-16"

// Flag
static int showVersion              = 0;     // just print version and quit?
static int showHelp                 = 0;     // just print help and quit?

// Argument constants for getopts
static const int ARG_log            = 256;
static const int ARG_error          = 257;
static const int ARG_readPairFile   = 258;
static const int ARG_outDGFile      = 259;
static const int ARG_minOverlap     = 260;
static const int ARG_multipleDG     = 261;
static const int ARG_inDGFile       = 262;
static const int ARG_outCollapseFile= 263;
static const int ARG_maxGap         = 264;
static const int ARG_maxDGOverhang  = 265;

// General Parameters
static string logFile = "";
static string errorFile = "";

// GenDuplexGroup Parameters
static string inReadPairFile = "";
static string outDGFile = "";
static int minOverlap = 5;
static bool multipleDG = false;

// CollapseDuplexGroup Parameters
static string inDGFile = "";
static string outCollapseFile = "";
static int maxGap = 10;
static int maxDGOverhang = 30;

/**
 * Print a detailed usage message to the provided output stream.
 */
static void printUsage(ostream& out) {
    out << "Usage: paris command [options] --log fileName --error fileName" << endl << endl
    << "    command:                    GenDuplexGroup / CollapseDuplexGroup" << endl
    << "        --log                   <String>"<< endl
    << "                                    File to Record the log Information" << endl
    << "        --error                 <String>"<< endl
    << "                                    File to Record the error Information" << endl
    << "        --help                  print this information and quit" << endl
    << "        --version               print version information and quit" << endl << endl
    << "    1. GenDuplexGroup options" << endl << endl
    << "        --inReadPairFile        <String>"<< endl
    << "                                    Sorted Read Pair File From Upstream PARIS pipeline" << endl
    << "        --outDGFile             <String>"<< endl
    << "                                    Output Duplex Group File" << endl
    << "        --minOverlap            <Interger>" << endl
    << "                                    minimum number of overlap nucleotides to define a read duplex (default: 5)" << endl
    << "        --multipleDG            <yes/no>" << endl
    << "                                    allow one read in multiple duplex groups (yes/no default: no)" << endl << endl
    << "    2. CollapseDuplexGroup options" << endl << endl
    << "        --inDGFile              <String>"<< endl
    << "                                    Sorted Duplex Group File From Upstream PARIS pipeline" << endl
    << "        --outCollapseFile       <String>"<< endl
    << "                                    Output Collapsed Duplex Group File" << endl
    << "        --maxGap                <Interger>" << endl
    << "                                    how many nucleotides allowed gapped during merging duplex group (default: 10)" << endl
    << "        --maxDGOverhang         <Interger>" << endl
    << "                                    how long the DG arm allowed (default: 30)" << endl
    ;
}

static const char *short_options = "hv";

/* https://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Options.html */
static struct option long_options[] = {
    {"help",            no_argument,            &showHelp,      'h'},
    {"version",         no_argument,            &showVersion,   'v'},
    
    {"log",             required_argument,      0,              ARG_log},
    {"error",           required_argument,      0,              ARG_error},
    
    {"inReadPairFile",  required_argument,      0,              ARG_readPairFile},
    {"outDGFile",       required_argument,      0,              ARG_outDGFile},
    {"minOverlap",      required_argument,      0,              ARG_minOverlap},
    {"multipleDG",      required_argument,      0,              ARG_multipleDG},
    
    {"inDGFile",         required_argument,     0,              ARG_inDGFile},
    {"outCollapseFile",  required_argument,     0,              ARG_outCollapseFile},
    {"maxGap",           required_argument,     0,              ARG_maxGap},
    {"maxDGOverhang",    required_argument,     0,              ARG_maxDGOverhang},
    {0, 0, 0, 0} // terminator
};

/**
 * Parse an int out of optarg and enforce that it be at least 'lower';
 * if it is less than 'lower', then output the given error message and
 * exit with an error and a usage message.
 */
template<typename T>
static int parseNumber(T lower, const char *errmsg) {
    char *endPtr= NULL;
    T t = (T)strtoll(optarg, &endPtr, 10);
    if (endPtr != NULL) {
        if (t < lower) {
            cerr << errmsg << endl;
            printUsage(cerr);
            exit(1);
        }
        return t;
    }
    cerr << errmsg << endl;
    printUsage(cerr);
    exit(1);
    return -1;
}

/**
 * Read command-line arguments
 */
static void parseOptions(int argc, char **argv) {
    int option_index = 0;
    int next_option;
    do {
        //cout << "Main " << endl;
        next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
       // cout << "optarg: " << next_option << "    " << optarg << endl;
        switch (next_option) {
            case ARG_log:
                logFile = optarg;
                break;
            case ARG_error:
                errorFile = optarg;
                break;
            
            case ARG_readPairFile:
                inReadPairFile = optarg;
                break;
            case ARG_outDGFile:
                outDGFile = optarg;
                break;
            case ARG_minOverlap:
                minOverlap = parseNumber<int>(1, "--minOverlap arg must be at least 1");
                break;
            case ARG_multipleDG:
                if(string(optarg) == "yes")
                    multipleDG = true;
                break;
            
            case ARG_inDGFile:
                inDGFile = optarg;
                break;
            case ARG_outCollapseFile:
                outCollapseFile = optarg;
                break;
            case ARG_maxGap:
                maxGap = parseNumber<int>(0, "--maxGap arg must be at least 0");
                break;
            case ARG_maxDGOverhang:
                maxDGOverhang = parseNumber<int>(10, "--maxDGOverhang arg must be at least 10");
                break;
        
            case -1: /* Done with options. */
                break;
            case 0:
                if (long_options[option_index].flag != 0)
                    break;	
            default: 
                cerr << "Unknown option: " << (char)next_option << endl;
                printUsage(cerr);
                exit(1);
        }
    } while(next_option != -1);
}



int main(int argc, char * argv[]) {
    
    
   // string infile;
   // vector<string> infiles;
   // string outfile;
    
    if (argc < 3) {
        printUsage(cerr);
        return -1;
    }
    
    //cout << string(argv[1]) << endl;
    if ( string(argv[1]) != "GenDuplexGroup" and string(argv[1]) != "CollapseDuplexGroup" )
    {
        cerr << "Specify Your Mode: GenDuplexGroup/CollapseDuplexGroup" << endl;
        return -1;
    }
    
  //  char **opt_argv = argv + 1;
  //  for(int i=0;i<argc-1;i++)
   //     cout << opt_argv[i] << endl;
    
   // return -1;
    
    parseOptions(argc-1, argv+1);
    string argv0 = argv[0];
    if(showVersion) {
        cout << "=================================" << endl;
        cout << "paris Version: " << PARIS_VERSION << endl;
        cout << "Built On: " << BUILD_HOST << endl;
        cout << "Built Time: " << BUILD_TIME << endl;
        cout << "=================================" << endl;
        //cout << "Source hash: " << EBWT_BUILD_HASH << endl;
        return 0;
    }
    
   // cout << "Main " << endl;
    
    if(showHelp) {
        printUsage(cerr);
        return 0;
    }
    
    string mode = argv[1];
    
    if (mode == "GenDuplexGroup")
    {
        //cout << "multipleDG: " << multipleDG << endl;
        if(inReadPairFile != "" and outDGFile != "")
            genDuplexGroup(inReadPairFile, outDGFile, logFile, errorFile, minOverlap, multipleDG);
        else
            cerr << "Error" << endl;
    }else if (mode == "CollapseDuplexGroup")
    {
        if(inDGFile != "" and outCollapseFile != "")
            collapseDuplexGroup(inDGFile, outCollapseFile, logFile, errorFile, maxGap, maxDGOverhang);
        else
            cerr << "Error" << endl;
    }
    
    return 0;
}
