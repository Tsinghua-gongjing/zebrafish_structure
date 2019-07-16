//
//  TimeFunc.hpp
//  GenDuplexGroup
//
//  Created by Lee on 2017/7/5.
//  Copyright © 2017年 Lee. All rights reserved.
//

#ifndef TimeFunc_hpp
#define TimeFunc_hpp


#include <stdio.h>
#include <sys/time.h>
#include <unistd.h>
#include <iostream>

using namespace std;

typedef long double(* microTransFunc)(long long);

long int getCurrentTime(); // get current microSeconds

long double microToMilli(long long microNumber);
long double microToSecond(long long microNumber);
long double microToMinute(long long microNumber);
long double microToHour(long long microNumber);
long double microToDay(long long microNumber);
long double microToWeek(long long microNumber);
long double microToMonth(long long microNumber);
long double microToYear(long long microNumber);

string MicroIntervalToHuman(long long microNumber);
string currentDate();

class ProcessAlert
{
public:
    ProcessAlert(long long total, long stepToAlert, int careRecentTimes);
    string alert(string other="");
    string finish();
    ~ProcessAlert();
private:
    //long long total;
    int careRecentTimes;
    long *recent;
    int recentIndex;
    long long total;
    long long current;
    long lastTime;
    long startTime;
    long stepToAlert;
    long long totalStep;
    
//flag
    bool full;
    
//private function
    long averageRecentTime();
};

#endif /* TimeFunc_hpp */


