#ifndef TIMESPAN_H_
#define TIMESPAN_H_

namespace
{
    static const long long TicksPerDay =  86400000000LL;
    static const long long TicksPerHour =  3600000000LL;
    static const long long TicksPerMinute =  60000000LL;
    static const long long TicksPerSecond =   1000000LL;
    static const long long TicksPerMillisecond = 1000LL;
    static const long long TicksPerMicrosecond =    1LL;

    static const long long UnixEpoch = 62135596800000000LL;

    static const long long MaxValueTicks = 315537897599999999LL;

    // 1582-Oct-15
    static const long long GregorianStart = 49916304000000000LL;
}

class TimeSpan
{
public:
    TimeSpan(long long ticks)
        : m_ticks(ticks)
    {
    }

    TimeSpan(int hours, int minutes, int seconds)
    {
        CalculateTicks(0, hours, minutes, seconds, 0);
    }

    TimeSpan(int days, int hours, int minutes, int seconds)
    {
        CalculateTicks(days, hours, minutes, seconds, 0);
    }

    TimeSpan(int days, int hours, int minutes, int seconds, int microseconds)
    {
        CalculateTicks(days, hours, minutes, seconds, microseconds);
    }

    TimeSpan Add(const TimeSpan& ts) const
    {
        return TimeSpan(m_ticks + ts.m_ticks);
    }
    
    TimeSpan Subtract(const TimeSpan& ts) const
    {
        return TimeSpan(m_ticks - ts.m_ticks);
    }

    int Compare(const TimeSpan& ts) const
    {
        int ret = 0;

        if (m_ticks < ts.m_ticks)
        {
            ret = -1;
        }
        if (m_ticks < ts.m_ticks)
        {
            ret = 1;
        }
        return ret;
    }

    bool Equals(const TimeSpan& ts) const
    {
        return m_ticks == ts.m_ticks;
    }

    int Days()
    {
        return m_ticks / TicksPerDay;
    }

    int Hours()
    {
        return (m_ticks % TicksPerDay / TicksPerHour);
    }

    int Minutes() const
    {
        return (m_ticks % TicksPerHour / TicksPerMinute);
    }

    int Seconds() const
    {
        return (m_ticks % TicksPerMinute / TicksPerSecond);
    }

    int Milliseconds() const
    {
        return (m_ticks % TicksPerSecond / TicksPerMillisecond);
    }
    
    int Microseconds() const
    {
        return (m_ticks % TicksPerSecond / TicksPerMicrosecond);
    }

    long long Ticks() const
    {
        return m_ticks;
    }

    double TotalDays() const
    {
        return m_ticks / static_cast<double>(TicksPerDay);
    }

    double TotalHours() const
    {
        return m_ticks / static_cast<double>(TicksPerHour);
    }

    double TotalMinutes() const
    {
        return m_ticks / static_cast<double>(TicksPerMinute);
    }

    double TotalSeconds() const
    {
        return m_ticks / static_cast<double>(TicksPerSecond);
    }
    
    double TotalMilliseconds() const
    {
        return m_ticks / static_cast<double>(TicksPerMillisecond);
    }
    
    double TotalMicroseconds() const
    {
        return m_ticks / static_cast<double>(TicksPerMicrosecond);
    }

private:
    long long m_ticks;

    void CalculateTicks(int days,
            int hours,
            int minutes,
            int seconds,
            int microseconds)
    {
        m_ticks = days * TicksPerDay +
            (hours * 3600LL + minutes * 60LL + seconds) * TicksPerSecond + 
            microseconds * TicksPerMicrosecond;
    }
};

inline TimeSpan operator+(const TimeSpan& ts1, const TimeSpan& ts2)
{
    return ts1.Add(ts2);
}

inline TimeSpan operator-(const TimeSpan& ts1, const TimeSpan& ts2)
{
    return ts1.Subtract(ts2);
}

inline bool operator==(const TimeSpan& ts1, TimeSpan& ts2)
{
    return ts1.Equals(ts2);
}

inline bool operator>(const TimeSpan& ts1, const TimeSpan& ts2)
{
    return (ts1.Compare(ts2) > 0);
}

inline bool operator>=(const TimeSpan& ts1, const TimeSpan& ts2)
{
    return (ts1.Compare(ts2) >= 0);
}

inline bool operator!=(const TimeSpan& ts1, const TimeSpan& ts2)
{
    return !ts1.Equals(ts2);
}

inline bool operator<(const TimeSpan& ts1, const TimeSpan& ts2)
{
    return (ts1.Compare(ts2) < 0);
}

inline bool operator<=(const TimeSpan& ts1, const TimeSpan& ts2)
{
    return (ts1.Compare(ts2) <= 0);
}

#endif
