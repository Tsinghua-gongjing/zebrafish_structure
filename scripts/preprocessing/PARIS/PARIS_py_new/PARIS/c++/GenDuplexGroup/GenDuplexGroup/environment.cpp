//
//  environment.cpp
//  GenDuplexGroup
//
//  Created by Lee on 2017/7/5.
//  Copyright © 2017年 Lee. All rights reserved.
//

#include "environment.hpp"



/*
 
 
 Split Function
 Time-consuming
 
 */

template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}
std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}



/*
 
 
 Read Class
 
 
 */

Read::Read(string chr1, int start1, int end1, char strand1,
           string chr2, int start2, int end2, char strand2,
           int readID)
{
    this->chr1 = chr1;
    this->start1 = start1;
    this->end1 = end1;
    this->strand1 = strand1;
    this->chr2 = chr2;
    this->start2 = start2;
    this->end2 = end2;
    this->strand2 = strand2;
    this->id = readID;
}

/*
 
 
 DuplexGroup Class
 
 
 */

DuplexGroup::DuplexGroup(const Read & read)
{
    this->chr1 = read.chr1;
    this->start1 = read.start1;
    this->end1 = read.end1;
    this->strand1 = read.strand1;
    this->chr2 = read.chr2;
    this->start2 = read.start2;
    this->end2 = read.end2;
    this->strand2 = read.strand2;
    this->support = 1;
    this->supportRead = to_string(read.id);
    this->collapsedTo = nullptr;
}

DuplexGroup::DuplexGroup(string chr1, int start1, int end1, char strand1, string chr2, int start2, int end2, char strand2, int support, string supportRead)
{
    this->chr1 = chr1;
    this->start1 = start1;
    this->end1 = end1;
    this->strand1 = strand1;
    this->chr2 = chr2;
    this->start2 = start2;
    this->end2 = end2;
    this->strand2 = strand2;
    this->support = support;
    this->supportRead = supportRead;
    this->collapsedTo = nullptr;
}

void DuplexGroup::addRead(const Read & read)
{
    this->start1 = max(this->start1, read.start1);
    this->end1 = min(this->end1, read.end1);
    
    this->start2 = max(this->start2, read.start2);
    this->end2 = min(this->end2, read.end2);
    
    this->support += 1;
    this->supportRead += ";"+to_string(read.id);
}

int DuplexGroup::overlapRead(const Read &read)
{
    // Suppose other_dg is after this one
    // -1: beyound the overlap
    // 0: no beyound but no overlap
    // other: has overlap
    if(chr1 != read.chr1 || chr2 != read.chr2 || strand1 != read.strand1 || strand2 != read.strand2)
        return -1;
    if(read.start1 > end1)
        return -1;
    if(read.start1 <= end1 and start2 <= read.end2 and read.start2 <= end2)
    {
        int overlap1 = min(end1, read.end1) - max(start1, read.start1);
        int overlap2 = min(end2, read.end2) - max(start2, read.start2);
        return min(overlap1, overlap2);
    }
    return 0;
}

int DuplexGroup::overlapDG(const DuplexGroup *dg, int maxGap=10, int maxTotal=30)
{
    int overlap = 0;
    if( end1 + maxGap < dg->start1 or chr1 != dg->chr1 or strand1 != dg->strand1 or chr2 != dg->chr2 or strand2 != dg->strand2 )
    {
        /*
        if( dg->supportRead == "23198;10015265;10015267;10015268" and supportRead=="10015266" )
        {
            cout << end1 << " " << dg->start2 << " " << int(end1 + maxGap < dg->start2) << endl;
        }
        */
         
        overlap = -1;
    }
    else{
        int gap2 = dg->start2 > start2 ? dg->start2 - end2 : start2 - dg->end2;
        int total1 = max(dg->end1, end1) - start1;
        int total2 = max(end2, dg->end2) - min(start2, dg->start2);
        if (total1 < maxTotal and total2 < maxTotal and gap2 < maxGap)
            overlap = 1;
    }
    return overlap;
}



void DuplexGroup::mergeDuplexGroup(DuplexGroup* dg)
{
    
    
    /*
    string combinedReads = supportRead + ";" + dg->supportRead;
    std::vector<std::string> reads = split(combinedReads, ';');
    
    string uniq_reads = "";
    PointerArray<string *> array;
    for(string read: reads)
    {
        string *current_read = new string(read);
        
        if(array.find(current_read) == -1)
        {
         //   cout << read << endl;
            array.append(current_read);
            uniq_reads += uniq_reads == "" ? read : ";" + read;
        }
    }
    
    */
    
    /*
    for(int i=0; i<reads.size(); i++)
    {
        
        if(array.find(current_read) == -1)
        {
            cout << reads[i] << endl;
            array.append(current_read);
            uniq_reads += uniq_reads == "" ? reads[i] : ";" + reads[i];
        }
    }
    */
    
    supportRead = supportRead + ";" + dg->supportRead;
    
   // supportRead = uniq_reads;
    
    start1 = min(start1, dg->start1);
    end1 = max(end1, dg->end1);
    
    start2 = min(start2, dg->start2);
    end2 = max(end2, dg->end2);
    
    support += dg->support;
}


















