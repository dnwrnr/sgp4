#ifndef JULIAN_H_
#define JULIAN_H_

#include "Globals.h"
#include "Timespan.h"

#include <ctime>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>

class Julian
{
public:
    Julian();

    Julian(const Julian& jul)
    {
        date_ = jul.date_;
    }

    Julian(const time_t t);

    /*
     * create julian date from year and day of year
     */
    Julian(int year, double day)
    {
        Initialize(year, day);
    }

    /*
     * create julian date from individual components
     * year: 2004
     * mon: 1-12
     * day: 1-31
     * hour: 0-23
     * min: 0-59
     * sec: 0-59.99
     */
    Julian(int year, int mon, int day, int hour, int min, double sec)
    {
        Initialize(year, mon, day, hour, min, sec);
    }

    ~Julian()
    {
    }

    // comparison operators
    bool operator==(const Julian &date) const;
    bool operator!=(const Julian &date) const;
    bool operator>(const Julian &date) const;
    bool operator<(const Julian &date) const;
    bool operator>=(const Julian &date) const;
    bool operator<=(const Julian &date) const;

    // assignment
    Julian& operator =(const Julian& b);
    Julian& operator =(const double b);
    // arithmetic
    Julian operator +(const Timespan& b) const;
    Julian operator -(const Timespan& b) const;
    Timespan operator -(const Julian& b) const;
    // compound assignment
    Julian & operator +=(const Timespan& b);
    Julian & operator -=(const Timespan& b);

    std::string ToString() const
    {
        std::stringstream ss;
        struct DateTimeComponents dt;
        ToGregorian(&dt);
        ss << std::right << std::fixed;
        ss << std::setprecision(6) << std::setfill('0');
        ss << std::setw(4) << dt.years << "-";
        ss << std::setw(2) << dt.months << "-";
        ss << std::setw(2) << dt.days << " ";
        ss << std::setw(2) << dt.hours << ":";
        ss << std::setw(2) << dt.minutes << ":";
        ss << std::setw(9) << dt.seconds << " UTC";
        return ss.str();
    }

    struct DateTimeComponents {

        DateTimeComponents()
        : years(0), months(0), days(0),
        hours(0), minutes(0), seconds(0.0) {
        }

        int years;
        int months;
        int days;
        int hours;
        int minutes;
        double seconds;
    };

    void ToGregorian(struct DateTimeComponents* datetime) const;

    time_t ToTime() const;
    double ToGreenwichSiderealTime() const;
    double ToLocalMeanSiderealTime(const double& lon) const;

    double FromJan1_00h_1900() const
    {
        return date_ - kEPOCH_JAN1_00H_1900;
    }

    double FromJan1_12h_1900() const
    {
        return date_ - kEPOCH_JAN1_12H_1900;
    }

    double FromJan1_12h_2000() const
    {
        return date_ - kEPOCH_JAN1_12H_2000;
    }

    double GetDate() const
    {
        return date_;
    }

    void AddDay(double day)
    {
        date_ += day;
    }

    void AddHour(double hr)
    {
        date_ += (hr / kHOURS_PER_DAY);
    }

    void AddMin(double min)
    {
        date_ += (min / kMINUTES_PER_DAY);
    }

    void AddSec(double sec) 
    {
        date_ += (sec / kSECONDS_PER_DAY);
    }

    static bool IsLeapYear(int y)
    {
        return (y % 4 == 0 && y % 100 != 0) || (y % 400 == 0);
    }

protected:
    void Initialize(int year, double day);
    void Initialize(int year, int mon, int day, int hour, int min, double sec);

    /*
     * the stored julian date
     */
    double date_;
};


inline std::ostream& operator<<(std::ostream& strm, const Julian& j)
{
    return strm << j.ToString();
}

#endif
