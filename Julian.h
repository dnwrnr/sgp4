#ifndef JULIAN_H_
#define JULIAN_H_

#include "Globals.h"

#include <ctime>

class Julian {
public:
    Julian();
    Julian(const Julian& jul);
    Julian(const time_t t);
    Julian(int year, double day);
    Julian(int year, int mon, int day, int hour, int min, double sec);

    /*
     * comparison operators
     */
    bool operator==(const Julian &date) const;
    bool operator!=(const Julian &date) const;
    bool operator>(const Julian &date) const;
    bool operator<(const Julian &date) const;
    bool operator>=(const Julian &date) const;
    bool operator<=(const Julian &date) const;

    ~Julian() {
    };

    time_t ToTime() const;
    double ToGreenwichSiderealTime() const;
    double ToLocalMeanSiderealTime(const double& lon) const;

    double FromJan1_00h_1900() const {
        return date_ - Globals::EPOCH_JAN1_00H_1900();
    }

    double FromJan1_12h_1900() const {
        return date_ - Globals::EPOCH_JAN1_12H_1900();
    }

    double FromJan1_12h_2000() const {
        return date_ - Globals::EPOCH_JAN1_12H_2000();
    }

    void GetComponent(int& year, int& month, double& dom) const;

    double GetDate() const {
        return date_;
    }

    void AddDay(double day) {
        date_ += day;
    }

    void AddHour(double hr) {
        date_ += (hr / Globals::HR_PER_DAY());
    }

    void AddMin(double min) {
        date_ += (min / Globals::MIN_PER_DAY());
    }

    void AddSec(double sec) {
        date_ += (sec / Globals::SEC_PER_DAY());
    }

    double SpanDay(const Julian& b) const {
        return date_ - b.date_;
    }

    double SpanHour(const Julian& b) const {
        return SpanDay(b) * Globals::HR_PER_DAY();
    }

    double SpanMin(const Julian& b) const {
        return SpanDay(b) * Globals::MIN_PER_DAY();
    }

    double SpanSec(const Julian& b) const {
        return SpanDay(b) * Globals::SEC_PER_DAY();
    }

    static bool IsLeapYear(int y) {
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

#endif
