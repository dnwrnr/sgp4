#include "Globals.h"
#include "Julian.h"
#include "Util.h"

#include <cmath>
#include <ctime>
#include <cassert>

#ifdef WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

Julian::Julian()
{
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

/*
 * create julian date given time_t value
 */
Julian::Julian(const time_t t)
{
    struct tm ptm;
#if WIN32
    if (gmtime_s(&ptm, &t))
    {
        assert(1);
    }
#else
    if (gmtime_r(&t, &ptm) == NULL)
    {
        assert(1);
    }
#endif
    int year = ptm.tm_year + 1900;
    double day = ptm.tm_yday + 1 +
            (ptm.tm_hour +
            ((ptm.tm_min +
            (ptm.tm_sec / 60.0)) / 60.0)) / 24.0;

    Initialize(year, day);
}

/*
 * comparison
 */
bool Julian::operator==(const Julian &date) const
{
    return date_ == date.date_ ? true : false;
}

bool Julian::operator!=(const Julian &date) const
{
    return !(*this == date);
}

bool Julian::operator>(const Julian &date) const
{
    return date_ > date.date_ ? true : false;
}

bool Julian::operator<(const Julian &date) const {

    return date_ < date.date_ ? true : false;
}

bool Julian::operator>=(const Julian &date) const
{
    return date_ >= date.date_ ? true : false;
}

bool Julian::operator<=(const Julian &date) const
{
    return date_ <= date.date_ ? true : false;
}

/*
 * assignment
 */
Julian& Julian::operator=(const Julian& b)
{
    if (this != &b) {
        date_ = b.date_;
    }
    return (*this);
}

Julian& Julian::operator=(const double b)
{
    date_ = b;
    return (*this);
}

/*
 * arithmetic
 */
Julian Julian::operator +(const Timespan& b) const
{
    return Julian(*this) += b;
}

Julian Julian::operator-(const Timespan& b) const
{
    return Julian(*this) -= b;
}

Timespan Julian::operator-(const Julian& b) const
{
    return Timespan(date_ - b.date_);
}

/*
 * compound assignment
 */
Julian & Julian::operator +=(const Timespan& b)
{
    date_ += b;
    return (*this);
}

Julian & Julian::operator -=(const Timespan& b)
{
    date_ -= b;
    return (*this);
}

/*
 * create julian date from year and day of year
 */
void Julian::Initialize(int year, double day)
{
    year--;

    int A = (year / 100);
    int B = 2 - A + (A / 4);

    double new_years = static_cast<int> (365.25 * year) +
            static_cast<int> (30.6001 * 14) +
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
void Julian::Initialize(int year, int mon, int day,
        int hour, int min, double sec)
{
    // Calculate N, the day of the year (1..366)
    int N;
    int F1 = (int) ((275.0 * mon) / 9.0);
    int F2 = (int) ((mon + 9.0) / 12.0);

    if (IsLeapYear(year))
    {
        // Leap year
        N = F1 - F2 + day - 30;
    }
    else
    {
        // Common year
        N = F1 - (2 * F2) + day - 30;
    }

    double dblDay = N + (hour + (min + (sec / 60.0)) / 60.0) / 24.0;

    Initialize(year, dblDay);
}

/*
 * converts time to time_t
 * note: resolution to seconds only
 */
time_t Julian::ToTime() const
{
    return static_cast<time_t> ((date_ - 2440587.5) * 86400.0);
}

/*
 * Greenwich Mean Sidereal Time
 */
double Julian::ToGreenwichSiderealTime() const
{
#if 0
    const double UT = fmod(jul + 0.5, 1.0);
    const double TU = (jul - 2451545.0 - UT) / 36525.0;

    double GMST = 24110.54841 + TU *
        (8640184.812866 + TU * (0.093104 - TU * 6.2e-06));

    GMST = fmod(GMST + SEC_PER_DAY * OMEGA_E * UT, SEC_PER_DAY);

    if (GMST < 0.0)
    {
        GMST += SEC_PER_DAY;  // "wrap" negative modulo value
    }

    return  (TWOPI * (GMST / SEC_PER_DAY));
#endif

    // tut1 = Julian centuries from 2000 Jan. 1 12h UT1
    // (since J2000 which is 2451545.0)
    // a Julian century is 36525 days
    const double tut1 = (date_ - 2451545.0) / 36525.0;
    const double tut2 = tut1 * tut1;
    const double tut3 = tut2 * tut1;

    // Rotation angle in arcseconds
    double theta = 67310.54841 + (876600.0 * 3600.0 + 8640184.812866) * tut1
        + 0.093104 * tut2 - 0.0000062 * tut3;

    // 360.0 / 86400.0 = 1.0 / 240.0
    theta = Util::WrapTwoPI(Util::DegreesToRadians(theta / 240.0));

    return theta;

#if 0
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
    const double c1p2p = C1 + kTWOPI;
    double gsto = Util::WrapTwoPI(THGR70 + C1 * ds70 + c1p2p * tfrac
            + ts70 * ts70 * FK5R);

    return gsto;
#endif
}

/*
 * Local Mean Sideral Time
 */
double Julian::ToLocalMeanSiderealTime(const double& lon) const
{
    return fmod(ToGreenwichSiderealTime() + lon, kTWOPI);
}

void Julian::ToGregorian(struct DateTimeComponents* datetime) const
{
    double jdAdj = GetDate() + 0.5;
    int Z = (int) jdAdj;
    double F = jdAdj - Z;

    int A = 0;

    if (Z < 2299161)
    {
        A = static_cast<int> (Z);
    }
    else
    {
        int a = static_cast<int> ((Z - 1867216.25) / 36524.25);
        A = static_cast<int> (Z + 1 + a - static_cast<int> (a / 4));
    }

    int B = A + 1524;
    int C = static_cast<int> ((B - 122.1) / 365.25);
    int D = static_cast<int> (365.25 * C);
    int E = static_cast<int> ((B - D) / 30.6001);

    datetime->hours = static_cast<int> (F * 24.0);
    F -= datetime->hours / 24.0;
    datetime->minutes = static_cast<int> (F * 1440.0);
    F -= datetime->minutes / 1440.0;
    datetime->seconds = F * 86400.0;

    datetime->days = B - D - static_cast<int> (30.6001 * E);
    datetime->months = E < 14 ? E - 1 : E - 13;
    datetime->years = datetime->months > 2 ? C - 4716 : C - 4715;
}
