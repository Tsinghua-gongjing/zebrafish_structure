//
//  DuplexGroup.hpp
//  GenDuplexGroup
//
//  Created by Lee on 2017/7/6.
//  Copyright © 2017年 Lee. All rights reserved.
//

#ifndef DuplexGroup_hpp
#define DuplexGroup_hpp


#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>

#include "TimeFunc.hpp"
#include "Files.hpp"
#include "environment.hpp"

using namespace std;

void genDuplexGroup(string read_pair_file, string duplex_group_file, string logFile, string errorFile, int OVERLAP=5, bool multiDG=false);
void collapseDuplexGroup(string sorted_duplex_group_file, string collapsed_dg_file, string logFile, string errorFile, int maxGap=10, int maxTotal=30);

#endif /* DuplexGroup_hpp */
