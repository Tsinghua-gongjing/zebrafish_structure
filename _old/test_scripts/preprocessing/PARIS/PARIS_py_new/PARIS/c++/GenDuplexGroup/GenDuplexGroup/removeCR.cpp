//
//  removeChimericReads.cpp
//  GenomeAnnotationResource
//
//  Created by Lee on 2017/7/17.
//
//


//#define __DEBUG

#define VERSION "1.0"
#define DATE "2017-7-18"
#define AUTHOR "Li Pan(hnsfyfyzlp@126.com)"

#define NORMAL 0
#define PARAM_ERROR 1
#define IO_ERROR 2

#define LONG long long

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <sstream>
using namespace std;

bool verbose = false;

LONG countFileLines(std::string fileName)
{
    int numLines = 0;
    ifstream in(fileName);
    std::string unused;
    while ( std::getline(in, unused) )
    ++numLines;
    in.close();
    return numLines;
}

int readJuncReads(istream &JUNC, vector<string> &reads, LONG fileLines)
{
    clog << "Start to read junction read names..." << endl;
    
    reads.clear();
    string line;
    LONG counting = 0;
    while( getline(JUNC, line) )
    {
#ifdef __DEBUG
        clog << "read line: " << line << endl;
#endif
        istringstream juncline(line);
        vector<string> juncLineItems;
        string item;
        while(juncline >> item)
            juncLineItems.push_back(item);

#ifdef __DEBUG
        clog << "juncLineItems[9]: " << juncLineItems[9] << endl;
#endif

        reads.push_back(juncLineItems[9]);
        
        ++counting;
        if( verbose && counting % 1000000 == 0)
        {
            char buffer[200];
            sprintf(buffer, "\treadJuncReads: have read %.2f%%", 100.0*counting/fileLines);
            //clog << "readJuncReads: have read " << static_cast<double>(counting)/fileLines << endl;
            clog << buffer << endl;
        }
    }
    
    clog << "Start to sort read names..." << endl;
    std::sort (reads.begin(), reads.end());
    
#ifdef __DEBUG
    
    for(decltype(reads.size()) i=0; i<10; i++)
        clog << "first " << i << "reads " << reads[i] << endl;
    
    for(auto i=reads.size()-10; i<reads.size(); i++)
        clog << "last " << i << "reads " << reads[i] << endl;
    
#endif
    
    return NORMAL;
}

int readUnmmappedAndOutput(istream &READS, const vector<string> &reads, ostream &OUT, LONG fileLines)
{
    
    clog << "Start to read and write unmmapped reads..." << endl;
    
    LONG idx = 0;
    string curLine; bool output = false;
    while( getline(READS, curLine) )
    {
        if(idx % 4 == 0)
        {
            curLine.erase(curLine.begin(), curLine.begin()+1);

#ifdef __DEBUG
            if(idx < 10)
                cout << "Cutted RNA ID: " << curLine;
#endif
            if(!std::binary_search(reads.begin(), reads.end(), curLine))
                output = true;
            else
                output = false;
            curLine = "@" + curLine;
        }
        if(output)
            OUT << curLine << endl;
        ++idx;
        if( verbose && idx % 1000000 == 0)
        {
            char buffer[200];
            sprintf(buffer, "\treadUnmmappedAndOutput: have read %.2f%%", 100.0*idx/fileLines);
            clog << buffer << endl;

            //clog << "readUnmmappedAndOutput: have read " << static_cast<double>(idx)/fileLines << endl;
        }
    }
    
    return NORMAL;
}

int main(int argc, char *argv[])
{
    if(argc != 4 && argc != 5)
    {
        cerr << "\tVersion: "<< VERSION << "\n\tDate: " << DATE << "\n\tAUTHOR: " << AUTHOR
        << "\nUsage: \n\tCommand Line: removeCR [-v] chimericJuncFile UnmmappedReadsFile outFilterReadsFile" << endl;
        cerr << "\t-v: verbose, output schedule.\n" << endl;
        return PARAM_ERROR;
    }
    
    int start_idx = 1;
    if(argc == 5)
    {
        if( string(argv[1]) != "-v" )
        {
            cerr << "Fatal Error: unrecognized option: " << argv[1] << endl;
            return PARAM_ERROR;
        }
        verbose = true;
        start_idx = 2;
    }
    
    string chimericJuncFile(argv[start_idx]);
    string UnmmappedReadsFile(argv[start_idx+1]);
    string outFilterReadsFile(argv[start_idx+2]);
    
    clog << "Input parameters..." << endl;
    clog << "\t ChimericJuncFile: " << chimericJuncFile << endl;
    clog << "\t UnmmappedReadsFile: " << UnmmappedReadsFile << endl;
    clog << "\t OutFilterReadsFile: " << outFilterReadsFile << endl;
    
    ifstream JUNC(chimericJuncFile, ifstream::in);
    ifstream READS(UnmmappedReadsFile, ifstream::in);
    ofstream OUT(outFilterReadsFile, ofstream::out);
    
    
    if(!JUNC)
    {
        cerr << "Fatal Error: open file " << chimericJuncFile << " failed" << endl;
        return IO_ERROR;
    }
    if(!READS)
    {
        cerr << "Fatal Error: open file " << UnmmappedReadsFile << " failed" << endl;
        return IO_ERROR;
    }
    if(!OUT)
    {
        cerr << "Fatal Error: open file " << outFilterReadsFile << " failed" << endl;
        return IO_ERROR;
    }
    
    clog << "start count file lines..." << endl;
    LONG junctionFileLines = countFileLines(chimericJuncFile);
    LONG unmmappedReadsFileLines = countFileLines(UnmmappedReadsFile);
    
    clog << "\t junctionFileLines: " << junctionFileLines << endl;
    clog << "\t unmmappedReadsFileLines: " << unmmappedReadsFileLines << endl;
    
    vector<string> read_names;
    readJuncReads(JUNC, read_names, junctionFileLines);
    readUnmmappedAndOutput(READS, read_names, OUT, unmmappedReadsFileLines);
    
    JUNC.close();
    READS.close();
    OUT.close();
    
    return NORMAL;
}












