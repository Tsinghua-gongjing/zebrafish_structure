//
//  Files.cpp
//  GenDuplexGroup
//
//  Created by Lee on 2017/7/5.
//  Copyright © 2017年 Lee. All rights reserved.
//

#include "Files.hpp"
#include <fstream>
using namespace std;


long long countFileLines(std::string fileName)
{
    int numLines = 0;
    ifstream in(fileName);
    std::string unused;
    while ( std::getline(in, unused) )
        ++numLines;
    in.close();
    return numLines;
}


