#include "Globals.h"
#include "Julian.h"

#include <cmath>
#include <ctime>
#include <cassert>
#include <cstring>

#ifdef WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

Julian::Julian() {
#ifdef WIN32
    SYSTEMTIME st;
    GetSystemTime(&st);
    Initialize(st.wYear,
            st.wMonth,
            st.wDay,
            st.wHour,
            st.wMinute,
            (double) st.wSecond + (double) st.wMilliseconds / 1000.0);
#else
    struct timeval tv;
    gettimeofday(&tv, NULL);
    struct tm gmt;
    gmtime_r(&tv.tv_sec, &gmt);
    Initialize(gmt.tm_year + 1900,
            gmt.tm_mon + 1,
            gmt.tm_mday,
            gmt.tm_hour,
            gmt.tm_min,
            (double) gmt.tm_sec + (double) tv.tv_usec / 1000000.0);
#endif
}

Julian::Julian(const Julian& jul) {
    date_ = jul.date_;
}

/*
 * create julian date given time_t value
 */
Julian::Julian(const time_t t) {
    struct tm ptm;
#if WIN32
    assert(gmtime_s(&ptm, &t));
#else
    if (gmtime_r(&t, &ptm) == NULL)
        assert(1);
#endif
    int year = ptm.tm_year + 1900;
    double day = ptm.tm_yday + 1 +
            (ptm.tm_hour +
            ((ptm.tm_min +
            (ptm.tm_sec / 60.0)) / 60.0)) / 24.0;

    Initialize(year, day);
}

/*
 * create julian date from year and day of year
 */
Julian::Julian(int year, double day) {
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
Julian::Julian(int year, int mon, int day, int hour, int min, double sec) {
    Initialize(year, mon, day, hour, min, sec);
}

/*
 * comparison
 */
bool Julian::operator==(const Julian &date) const {
    return date_ == date.date_ ? true : false;
}

bool Julian::operator!=(const Julian &date) const {
    return date_ == date.date_ ? false : true;
}

bool Julian::operator>(const Julian &date) const {
    return date_ > date.date_ ? true : false;
}

bool Julian::operator<(const Julian &date) const {
    return date_ < date.date_ ? true : false;
}

bool Julian::operator>=(const Julian &date) const {
    return date_ >= date.date_ ? true : false;
}

bool Julian::operator<=(const Julian &date) const {
    return date_ <= date.date_ ? true : false;
}

/*
 * assignment
 */
Julian& Julian::operator=(const Julian& b) {

    if (this != &b) {
        date_ = b.date_;
    }
    return (*this);
}

Julian& Julian::operator=(const double b) {

    date_ = b;
    return (*this);
}

/*
 * arithmetic
 */
Julian Julian::operator+(const Timespan& b) const {
    Julian result(*this);
    result.date_ += b.GetTotalDays();
    return result;
}

Julian Julian::operator-(const Timespan& b) const {
    Julian result(*this);
    result.date_ -= b.GetTotalDays();
    return result;
}

/*
 * create julian date from year and day of year
 */
void Julian::Initialize(int year, double day) {

    // Now calculate Julian date
    year--;

    // Centuries are not leap years unless they divide by 400
    int A = (year / 100);
    int B = 2 - A + (A / 4);

    // 1720994.5 = Oct 30, year -1
    double new_years = (int) (365.25 * year) +
            (int) (30.6001 * 14) +
            1720994.5 + B;

    date_ = new_years + day;
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
void Julian::Initialize(int year, int mon, int day, int hour, int min, double sec) {
    // Calculate N, the day of the year (1..366)
    int N;
    int F1 = (int) ((275.0 * mon) / 9.0);
    int F2 = (int) ((mon + 9.0) / 12.0);

    if (IsLeapYear(year)) {
        // Leap year
        N = F1 - F2 + day - 30;
    } else {
        // Common year
        N = F1 - (2 * F2) + day - 30;
    }

    double dblDay = N + (hour + (min + (sec / 60.0)) / 60.0) / 24.0;

    Initialize(year, dblDay);
}

/*
 * gets year, month and day of month from julian date
 */
void Julian::GetComponent(int& year, int& month, double& dom) const {

    double jdAdj = GetDate() + 0.5;
    int Z = (int) jdAdj; // integer part
    double F = jdAdj - Z; // fractional part
    double alpha = (int) ((Z - 1867216.25) / 36524.25);
    double A = Z + 1 + alpha - (int) (alpha / 4.0);
    double B = A + 1524.0;
    int C = (int) ((B - 122.1) / 365.25);
    int D = (int) (C * 365.25);
    int E = (int) ((B - D) / 30.6001);

    dom = B - D - (int) (E * 30.6001) + F;
    month = (E < 13.5) ? (E - 1) : (E - 13);
    year = (month > 2.5) ? (C - 4716) : (C - 4715);
}

/*
 * converts time to time_t
 * note: resolution to seconds only
 */
time_t Julian::ToTime() const {

    int nYear;
    int nMonth;
    double dblDay;

    GetComponent(nYear, nMonth, dblDay);

    // dblDay is the fractional Julian Day (i.e., 29.5577).
    // Save the whole number day in nDOM and convert dblDay to
    // the fractional portion of day.
    int nDOM = (int) dblDay;

    dblDay = dblDay - nDOM;

    const int SEC_PER_MIN = 60;
    const int SEC_PER_HR = 60 * SEC_PER_MIN;
    const int SEC_PER_DAY = 24 * SEC_PER_HR;

    int secs = (int) ((dblDay * SEC_PER_DAY) + 0.5);

    // Create a "struct tm" type.
    // NOTE:
    // The "struct tm" type has a 1-second resolution. Any fractional
    // component of the "seconds" time value is discarded.
    struct tm tGMT;
    memset(&tGMT, 0, sizeof (tGMT));

    tGMT.tm_year = nYear - 1900; // 2001 is 101
    tGMT.tm_mon = nMonth - 1; // January is 0
    tGMT.tm_mday = nDOM; // First day is 1
    tGMT.tm_hour = secs / SEC_PER_HR;
    tGMT.tm_min = (secs % SEC_PER_HR) / SEC_PER_MIN;
    tGMT.tm_sec = (secs % SEC_PER_HR) % SEC_PER_MIN;
    tGMT.tm_isdst = 0; // No conversion desired

    time_t tEpoch = mktime(&tGMT);

    if (tEpoch != -1) {
        // Valid time_t value returned from mktime().
        // mktime() expects a local time which means that tEpoch now needs
        // to be adjusted by the difference between this time zone and GMT.
        //		tEpoch = tEpoch - _timezone;
    }

    return tEpoch;
}

/*
 * Greenwich Mean Sidereal Time
 */
double Julian::ToGreenwichSiderealTime() const {
#if 0
    double theta;
    double tut1;

    // tut1 = Julian centuries from 2000 Jan. 1 12h UT1 (since J2000 which is 2451545.0)
    // a Julian century is 36525 days
    tut1 = (date_ - 2451545.0) / 36525.0;

    // Rotation angle in arcseconds
    theta = 67310.54841 + (876600.0 * 3600.0 + 8640184.812866) * tut1 + 0.093104 * pow(tut1, 2.0) - 0.0000062 * pow(tut1, 3.0);

    // 360.0 / 86400.0 = 1.0 / 240.0
    theta = fmod(Globals::Deg2Rad(theta / 240.0), Globals::TWOPI());

    /*
     * check quadrants
     */
    if (theta < 0.0)
        theta += Globals::TWOPI();

    return theta;
#endif

    static const double C1 = 1.72027916940703639e-2;
    static const double THGR70 = 1.7321343856509374;
    static const double FK5R = 5.07551419432269442e-15;

    /*
     * get integer number of days from 0 jan 1970
     */
    const double ts70 = date_ - 2433281.5 - 7305.0;
    const double ds70 = floor(ts70 + 1.0e-8);
    const double tfrac = ts70 - ds70;
    /*
     * find greenwich location at epoch
     */
    const double c1p2p = C1 + Globals::TWOPI();
    double gsto = fmod(THGR70 + C1 * ds70 + c1p2p * tfrac + ts70 * ts70 * FK5R, Globals::TWOPI());
    if (gsto < 0.0)
        gsto = gsto + Globals::TWOPI();

    return gsto;
}

/*
 * Local Mean Sideral Time
 */
double Julian::ToLocalMeanSiderealTime(const double& lon) const {
    return fmod(ToGreenwichSiderealTime() + lon, Globals::TWOPI());
}

void Julian::GetDateTime(struct DateTimeComponents* datetime) const {
    
    double jdAdj = GetDate() + 0.5;
    int Z = (int) jdAdj;
    double F = jdAdj - Z;

    int A = 0;

    if (Z < 2299161) {
        A = static_cast<int>(Z);
    } else {
        int a = static_cast<int>((Z - 1867216.25) / 36524.25);
        A = static_cast<int>(Z + 1 + a - static_cast<int>(a / 4));
    }

    int B = A + 1524;
    int C = static_cast<int>((B - 122.1) / 365.25);
    int D = static_cast<int>(365.25 * C);
    int E = static_cast<int>((B - D) / 30.6001);

    datetime->hours = static_cast<int>(F * 24.0);
    F -= datetime->hours / 24.0;
    datetime->minutes = static_cast<int>(F * 1440.0);
    F -= datetime->minutes / 1440.0;
    datetime->seconds = F * 86400.0;

    datetime->days = B - D - static_cast<int>(30.6001 * E);
    datetime->months = E < 14 ? E - 1 : E - 13;
    datetime->years = datetime->months > 2 ? C - 4716 : C - 4715;
 }
