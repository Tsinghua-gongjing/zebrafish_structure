//
//  DuplexGroup.cpp
//  GenDuplexGroup
//
//  Created by Lee on 2017/7/6.
//  Copyright © 2017年 Lee. All rights reserved.
//

#include "DuplexGroup.hpp"


void genDuplexGroup(string read_pair_file, string duplex_group_file, string logFile, string errorFile, int OVERLAP, bool multiDG)
{
    ofstream LOG(logFile, ios_base::app);
    ofstream ERR(errorFile, ios_base::app);
    
    LOG << "Start to genDuplexGroup " << read_pair_file << "  ===>  " << duplex_group_file << endl;
    
    long long fileLines = countFileLines(read_pair_file);
    ProcessAlert alert(fileLines, 100000, 10);
    
    ifstream READPAIR(read_pair_file);
    char strand1, strand2;
    string chr1, chr2;
    int start1, end1, id1, score1, start2, end2, id2, score2;
    
    int counter = 0;
    int firstPossible = 0;
    int dgCount = 0;
    PointerArray<DuplexGroup *> dgArray;
    while(READPAIR>>chr1>>start1>>end1>>id1>>score1>>strand1>>chr2>>start2>>end2>>id2>>score2>>strand2)
    {
        char buffer[200];
        sprintf(buffer, "firstPossible: %d, duplexCount: %d", firstPossible, dgCount);
        LOG << alert.alert(string(buffer));
        
        Read read(chr1, start1, end1, strand1, chr2, start2, end2, strand2, id1);
        
        int nonOverlapped = 1;
        int lastDGoverlapped = 0;
        for(int i=firstPossible;i<dgCount;i++)
        {
            int overlap = dgArray[i]->overlapRead(read);
            //cout << overlap << endl;
            if(overlap >= OVERLAP)
            {
                lastDGoverlapped = 1;
                nonOverlapped = 0;
                dgArray[i]->addRead(read);
                //cout << "find one" << endl;
                if (! multiDG) break;
            }else if (overlap == -1){
                if(not lastDGoverlapped)
                    firstPossible = i + 1;
            }else
            {
                lastDGoverlapped = 1;
            }
        }
        if (nonOverlapped)
        {
            DuplexGroup *dg = new DuplexGroup(read);
            dgArray.append( dg );
            dgCount++;
        }
        counter++;
    }
    LOG << alert.finish();
    READPAIR.close();
    
    ofstream DUPLEXGROUP(duplex_group_file);
    for(int i=0; i<dgArray.arrayLen(); i++)
        DUPLEXGROUP << dgArray[i]->chr1 << "\t" << dgArray[i]->start1 << "\t" << dgArray[i]->end1 << "\t" << dgArray[i]->strand1 << "\t" << dgArray[i]->chr2 << "\t" << dgArray[i]->start2 << "\t" << dgArray[i]->end2 << "\t" << dgArray[i]->strand2 << "\t" << dgArray[i]->supportRead << "\t" << dgArray[i]->support << endl;
    DUPLEXGROUP.close();
    
    LOG.close();
    ERR.close();
}

void collapseDuplexGroup(string sorted_duplex_group_file, string collapsed_dg_file, string logFile, string errorFile, int maxGap, int maxTotal)
{
    ofstream LOG(logFile, ios_base::app);
    ofstream ERR(errorFile, ios_base::app);
    
    LOG << "Start to collapseDuplexGroup " << sorted_duplex_group_file << "  ===>  " << collapsed_dg_file << endl;
    
    int merged_dgCount = 0;
    int firstPossible=0;
    
    PointerArray<DuplexGroup *> dgArray;
    long long fileLines = countFileLines(sorted_duplex_group_file);
    ProcessAlert alert(fileLines, 100000, 10);
    
    ifstream DUPLEXGROUP(sorted_duplex_group_file);
    char strand1, strand2;
    string chr1, chr2, supportReads;
    int start1, end1, start2, end2, support;
    
    while(DUPLEXGROUP>>chr1>>start1>>end1>>strand1>>chr2>>start2>>end2>>strand2>>supportReads>>support)
    {
        char buffer[200];
        sprintf(buffer, "firstPossible: %d, merged_dgCount: %d", firstPossible, merged_dgCount);
        LOG << alert.alert(string(buffer));
        
        DuplexGroup *dg = new DuplexGroup(chr1, start1, end1, strand1, chr2, start2,end2, strand2, support, supportReads);
        //dgArray.append(dg);
        
        //if(supportReads == "23198;10015265;10015267;10015268")
        //  cout << dgArray[dgArray.arrayLen()-1]->supportRead << endl;
        
        int lastDGovarlapped = 0; bool merged = false;
        for(int idx=firstPossible; idx<dgArray.arrayLen();idx++)
        {
            
            int overlapped = dgArray[idx]->overlapDG(dg, maxGap, maxTotal);
            /*
             if(supportReads == "23198;10015265;10015267;10015268" and idx==dgArray.arrayLen()-1)
             {
             cout << dgArray[dgArray.arrayLen()-1]->supportRead << endl;
             cout <<overlapped << endl;
             }
             */
            if (overlapped == -1){
                if(not lastDGovarlapped)
                    firstPossible = idx + 1;
            }
            else if(overlapped > 0){
                lastDGovarlapped = 1;
                dgArray[idx]->mergeDuplexGroup(dg);
                //dgArray.del(dgArray.arrayLen()-1);
                merged = true;
                delete dg;
                merged_dgCount++;
                break;
            }
            else{
                lastDGovarlapped = 1;
            }
        }
        if(not merged)
            dgArray.append(dg);
        //cout << merged_dgCount << endl;
    }
    LOG << alert.finish();
    DUPLEXGROUP.close();
    
    ofstream COLLAPSE(collapsed_dg_file);
    for(int i=0; i<dgArray.arrayLen(); i++)
        COLLAPSE << dgArray[i]->chr1 << "\t" << dgArray[i]->start1 << "\t" << dgArray[i]->end1 << "\t" << dgArray[i]->strand1 << "\t" << dgArray[i]->chr2 << "\t" << dgArray[i]->start2 << "\t" << dgArray[i]->end2 << "\t" << dgArray[i]->strand2 << "\t" << dgArray[i]->supportRead << "\t" << dgArray[i]->support << "\n";
    COLLAPSE.close();
    
    LOG.close();
    ERR.close();
}


