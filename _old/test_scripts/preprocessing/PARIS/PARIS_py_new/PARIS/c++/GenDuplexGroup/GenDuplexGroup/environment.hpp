//
//  environment.hpp
//  GenDuplexGroup
//
//  Created by Lee on 2017/7/5.
//  Copyright © 2017年 Lee. All rights reserved.
//

#ifndef environment_hpp
#define environment_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include <sstream>
using namespace std;

std::vector<std::string> split(const std::string &s, char delim);


class Read
{
public:
    Read(string chr1, int start1, int end1, char strand1,
         string chr2, int start2, int end2, char strand2,
         int readID);
    
public:
    string chr1, chr2;
    char strand1, strand2;
    int start1, end1, start2, end2;
    int id;
};



class DuplexGroup
{
public:
    DuplexGroup(const Read &);
    DuplexGroup(string chr1, int start1, int end1, char strand1, string chr2, int start2, int end2, char strand2, int support, string supportRead);
    void addRead(const Read &);
    int overlapRead(const Read &);
    int overlapDG(const DuplexGroup*, int, int);
    void mergeDuplexGroup(DuplexGroup*);
    //DuplexGroup *getCollapsedDG();
    
public:
    string chr1, chr2;
    char strand1, strand2;
    int start1, end1, start2, end2;
    int support;
    string supportRead;
    DuplexGroup *collapsedTo;
};




template <class T>
class PointerArray
{
public:
    PointerArray()
    {
        array = (T*)malloc(100*sizeof(T));
        content = 100;
        lastElem = 0;
    }
    ~PointerArray()
    {
        for(int i=0;i<lastElem;i++)
            delete array[i];
        delete array;
    }
    void append(T node)
    {
        /*
        if (lastElem + 1 == content)
        {
            array = (T*)realloc(array, (content+100)*sizeof(T ));
            content += 100;
        }
        array[lastElem] = node;
        lastElem += 1;
        */
        insert(lastElem, node);
    }
    void del(int elem)
    {
        if(elem < lastElem)
        {
            delete array[elem];
            for(int i=elem; i<lastElem-1; i++)
                array[i] = array[i+1];
            lastElem--;
        }
    }
    void insert(int elem, T node)
    {
        if(elem > lastElem)
            return;
        if (lastElem + 1 == content)
        {
            array = (T*)realloc(array, (content+100)*sizeof(T ));
            content += 100;
        }
        
        if(elem == lastElem)
        {
            array[lastElem] = node;
            lastElem ++;
            return;
        }
        
        for (int i=lastElem-1; i>=elem; i--)
            array[i+1] = array[i];
        
        array[elem] = node;
        lastElem++;
    }
    void clearAll()
    {
        for(int i=0; i<lastElem; i++)
            delete array[i];
        lastElem = 0;
    }
    
    int find(T obj, int elem=0)
    {
        for(int i=elem; i<lastElem; i++)
            if(*obj == *array[i])
                return i;
        return -1;
    }
    
    T operator[](int elem){
        if(elem < lastElem)
            return array[elem];
        else
            return nullptr;
    }
    int arrayLen(){ return lastElem; }

public:
    T *array;
    int content;
    int lastElem;
};



/* string as T is unstable */
template <class T>
class Array
{
public:
    Array()
    {
        array = (T*)malloc(100*sizeof(T));
        content = 100;
        lastElem = 0;
    }
    ~Array()
    {
        delete array;
    }
    void append(T node)
    {
        /*
        if (lastElem + 1 == content)
        {
            array = (T*)realloc(array, (content+100)*sizeof(T ));
            content += 100;
        }
        array[lastElem] = node;
        lastElem += 1;
        */
        insert(lastElem, node);
    }
    void del(int elem)
    {
        if(elem < lastElem)
        {
            //delete array[elem];
            for(int i=elem; i<lastElem-1; i++)
                array[i] = array[i+1];
            lastElem--;
        }
    }
    
    void insert(int elem, T node)
    {
        if(elem > lastElem)
            return;
        
        if (lastElem + 1 == content)
        {
            array = (T*)realloc(array, (content+100)*sizeof(T));
            content += 100;
        }
        
        //cout << elem << " " << lastElem << endl;
        if(elem == lastElem)
        {
            //cout << node << endl;
            array[lastElem] = node;
            lastElem ++;
            return;
            //cout << " out " << endl;
        }

        for (int i=lastElem-1; i>=elem; i--)
            array[i+1] = array[i];

        //cout << "Why?" << endl;
        array[elem] = node;
        lastElem++;
    }
    
    void clearAll()
    {
        lastElem = 0;
    }
    
    int find(T obj, int elem=0)
    {
        for(int i=elem; i<lastElem; i++)
            if(obj == array[i])
                return i;
        return -1;
    }
    
    T operator[](int elem){
        if(elem < lastElem)
            return array[elem];
        else
            return 0;
    }
    
    int arrayLen(){ return lastElem; }
    
public:
    T *array;
    int content;
    int lastElem;
};



#endif /* environment_hpp */
