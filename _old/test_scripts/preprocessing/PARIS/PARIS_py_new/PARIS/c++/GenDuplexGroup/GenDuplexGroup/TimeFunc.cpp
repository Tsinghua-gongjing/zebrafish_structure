//
//
//  GenDuplexGroup
//
//  Created by Lee on 2017/7/5.
//  Copyright © 2017年 Lee. All rights reserved.
//

#include "TimeFunc.hpp"
#include <sstream>

long int getCurrentTime()
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    long int ms = tp.tv_sec * 1000000 + tp.tv_usec;
    return ms;
}

long double microToMilli(long long microNumber)
{
    return double(microNumber) / 1000;
}

long double microToSecond(long long microNumber)
{
    return double(microNumber) / 1000 / 1000;
}

long double microToMinute(long long microNumber)
{
    return double(microNumber) / 1000 / 1000 / 60;
}

long double microToHour(long long microNumber)
{
    return double(microNumber) / 1000 / 1000 / 60 / 60;
}

long double microToDay(long long microNumber)
{
    return double(microNumber) / 1000 / 1000 / 60 / 60 / 24;
}

long double microToWeek(long long microNumber)
{
    return double(microNumber) / 1000 / 1000 / 60 / 60 / 24 / 7;
}

long double microToMonth(long long microNumber)
{
    return double(microNumber) / 1000 / 1000 / 60 / 60 / 24 / 7 / 4;
}

long double microToYear(long long microNumber)
{
    return double(microNumber) / 1000 / 1000 / 60 / 60 / 24 / 7 / 4 / 12;
}

long double microConst(long long microNumber)
{
    return double(microNumber);
}


/*
 cout << MicroIntervalToHuman(1LL) << endl;   //1.000 microSeconds
 cout << MicroIntervalToHuman(1000LL) << endl; //1.000 milliSeconds
 cout << MicroIntervalToHuman(1000000LL) << endl; //1.000 Seconds
 cout << MicroIntervalToHuman(60000000LL) << endl; //1.000 Minutes
 cout << MicroIntervalToHuman(60000000LL * 60) << endl; //1.000 Hours
 cout << MicroIntervalToHuman(60000000LL * 60 * 24) << endl; //1.000 Days
 cout << MicroIntervalToHuman(60000000LL * 60 * 24 * 7) << endl; //1.000 Weeks
 cout << MicroIntervalToHuman(60000000LL * 60 * 24 * 7 * 4) << endl; //1.000 Months
 cout << MicroIntervalToHuman(60000000LL * 60 * 24 * 7 * 4 * 12) << endl; //1.000 Years
 cout << MicroIntervalToHuman(60000000LL * 60 * 24 * 7 * 4 * 12 * 12) << endl; //12.000 Years
 */

string MicroIntervalToHuman(long long microNumber)
{
    if (microNumber <= 0) {
        return "0 microSeconds";
    }
    microTransFunc FunctionList[9] = {microConst, microToMilli, microToSecond, microToMinute, microToHour, microToDay, microToWeek, microToMonth, microToYear};
    string postFix[9] = {"microSeconds", "milliSeconds", "Seconds", "Minutes", "Hours", "Days", "Weeks", "Months", "Years"};
    
    int idx = 0;
    while( idx < 9 and FunctionList[idx](microNumber) >= 1.0 )
    {
        //cout << idx << "  " << FunctionList[idx](microNumber) << endl;
        idx++;
    }
    //idx = idx-1<0 ? 0 : idx-1;
    double humanTime = FunctionList[idx-1](microNumber);
    char buffer[50];
    sprintf(buffer, "%.3f", humanTime);
    return string(buffer) + " " + postFix[idx-1];
}

string currentDate()
{
    time_t timep;
    struct tm *p;
    time(&timep);
    p = localtime(&timep);
    
    char buffer[200];
    sprintf(buffer, "Current Time: %d-%d-%d; %d:%d:%d", 1900+p->tm_year, 1+p->tm_mon, p->tm_mday, p->tm_hour, p->tm_min, p->tm_sec);
    return string(buffer);
}

/*
 
 class ProcessAlert
 {
 public:
     ProcessAlert(long long total, long stepToAlert, int careRecentTimes);
     string alert(string other="");
     string finish();
     ~ProcessAlert();
 private:
     long long total;
     int careRecentTimes;
     long *recent;
    int current;
     long lastTime;
     long startTime;
     long stepToAlert;
     long long totalStep;
     int recentIndex;
 };
 
 */

ProcessAlert::ProcessAlert(long long total, long stepToAlert, int careRecentTimes)
{
    //flag
    full = false;
    
    //useful variables
    this->total = total;
    this->stepToAlert = stepToAlert;
    this->careRecentTimes = careRecentTimes;
    this->recent = new long[careRecentTimes];
    for(int i=0;i<careRecentTimes;i++)
        this->recent[i] = 0.0;
    this->lastTime = this->startTime = getCurrentTime();
    this->current = 0;
    this->totalStep = total / stepToAlert;
    this->recentIndex = 0;
}

long ProcessAlert::averageRecentTime()
{
    if ( not full)
    {
        // round to 0 when full
        if ( recentIndex != 0 )
            return -1;
        else
            full = true;
    }
    long sum = 0;
    for (int i=0;i<careRecentTimes;i++)
        sum += recent[i];
    return sum/careRecentTimes;
}

string ProcessAlert::alert(string other)
{
    if(current > total)
    {
        return "Error: Out of limitation\n";
    }
    current++;
    if ( current % stepToAlert == 0 )
    {
        double runnedPercent = 1.0 * double(current) / total;
        int leftStep = (1 - runnedPercent) * totalStep;
        long nowTime = getCurrentTime();
        long recentInterval = nowTime - lastTime;
        long speed = 1000000L * stepToAlert / (recentInterval+1);
        lastTime = nowTime;
        //cout << "recentIndex" << recentIndex << endl;
        recent[recentIndex] = recentInterval;
        recentIndex = ++recentIndex % careRecentTimes;
        long aveRecentTime = averageRecentTime();
      //  cout << "aveRecentTime: " << aveRecentTime << endl;
      //  cout << "left step: " << leftStep << endl;
        
        char buffer[200];
        if(aveRecentTime != -1)
        {
            long long leftTime = aveRecentTime * leftStep;
            string humanLeftTime = MicroIntervalToHuman(leftTime);
            sprintf(buffer, "Have read %lld lines(%.2f%%)...estimate %s left. %s\n\t"
                    "Speed: %ld lines/second\n\t%s\n",current, 100.0*current/total, humanLeftTime.c_str(), other.c_str(), speed, currentDate().c_str());
           // for(int i=0; i<careRecentTimes; i++)
           //     cout << recent[i] << "\t";
           // cout << endl;
        }
        else{
            sprintf(buffer, "Have read %lld lines(%.2f%%)... %s\n\t"
                    "Speed: %ld lines/second\n\t%s\n",current, 100.0*current/total, other.c_str(), speed, currentDate().c_str());
        }
        return string(buffer);
    }
    return "";
}


string ProcessAlert::finish()
{
    long nowTime = getCurrentTime();
    string total_time = MicroIntervalToHuman(nowTime - startTime);
    long ave_speed = 1000000 * current / (nowTime - startTime+1);
    char buffer[200];
    sprintf(buffer, "Task finished seccessfully: total time: %s; average speed: %ld lines/second.\n\t%s\n", total_time.c_str(), ave_speed, currentDate().c_str());
    return string(buffer);
}

ProcessAlert::~ProcessAlert()
{
    delete recent;
}
















