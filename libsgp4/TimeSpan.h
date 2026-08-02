/*
 * Copyright 2013 Daniel Warner <contact@danrw.com>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#pragma once

#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace libsgp4
{

namespace
{
    static const int64_t kTICKS_PER_DAY = 86400000000LL;
    static const int64_t kTICKS_PER_HOUR = 3600000000LL;
    static const int64_t kTICKS_PER_MINUTE = 60000000LL;
    static const int64_t kTICKS_PER_SECOND = 1000000LL;
    static const int64_t kTICKS_PER_MILLISECOND = 1000LL;
    static const int64_t kTICKS_PER_MICROSECOND = 1LL;

    static const int64_t kUNIX_EPOCH = 62135596800000000LL;

    static const int64_t kMAX_VALUE_TICKS = 315537897599999999LL;

    // 1582-Oct-15
    static const int64_t kGREGORIAN_START = 49916304000000000LL;
} // namespace

/**
 * @brief Represents a time interval.
 *
 * Represents a time interval (duration/elapsed) that is measured as a positive
 * or negative number of days, hours, minutes, seconds, and fractions
 * of a second.
 */
class TimeSpan
{
public:
    explicit TimeSpan(int64_t ticks)
        : mTicks(ticks)
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
        return TimeSpan(mTicks + ts.mTicks);
    }

    TimeSpan Subtract(const TimeSpan& ts) const
    {
        return TimeSpan(mTicks - ts.mTicks);
    }

    int Compare(const TimeSpan& ts) const
    {
        int ret = 0;

        if (mTicks < ts.mTicks)
        {
            ret = -1;
        }
        if (mTicks > ts.mTicks)
        {
            ret = 1;
        }
        return ret;
    }

    bool Equals(const TimeSpan& ts) const
    {
        return mTicks == ts.mTicks;
    }

    int Days() const
    {
        return static_cast<int>(mTicks / kTICKS_PER_DAY);
    }

    int Hours() const
    {
        return static_cast<int>(mTicks % kTICKS_PER_DAY / kTICKS_PER_HOUR);
    }

    int Minutes() const
    {
        return static_cast<int>(mTicks % kTICKS_PER_HOUR / kTICKS_PER_MINUTE);
    }

    int Seconds() const
    {
        return static_cast<int>(mTicks % kTICKS_PER_MINUTE / kTICKS_PER_SECOND);
    }

    int Milliseconds() const
    {
        return static_cast<int>(mTicks % kTICKS_PER_SECOND / kTICKS_PER_MILLISECOND);
    }

    int Microseconds() const
    {
        return static_cast<int>(mTicks % kTICKS_PER_SECOND / kTICKS_PER_MICROSECOND);
    }

    int64_t Ticks() const
    {
        return mTicks;
    }

    double TotalDays() const
    {
        return static_cast<double>(mTicks) / kTICKS_PER_DAY;
    }

    double TotalHours() const
    {
        return static_cast<double>(mTicks) / kTICKS_PER_HOUR;
    }

    double TotalMinutes() const
    {
        return static_cast<double>(mTicks) / kTICKS_PER_MINUTE;
    }

    double TotalSeconds() const
    {
        return static_cast<double>(mTicks) / kTICKS_PER_SECOND;
    }

    double TotalMilliseconds() const
    {
        return static_cast<double>(mTicks) / kTICKS_PER_MILLISECOND;
    }

    double TotalMicroseconds() const
    {
        return static_cast<double>(mTicks) / kTICKS_PER_MICROSECOND;
    }

    std::string ToString() const
    {
        std::stringstream ss;

        ss << std::right << std::setfill('0');

        if (mTicks < 0)
        {
            ss << '-';
        }

        if (Days() != 0)
        {
            ss << std::setw(2) << std::abs(Days()) << '.';
        }

        ss << std::setw(2) << std::abs(Hours()) << ':';
        ss << std::setw(2) << std::abs(Minutes()) << ':';
        ss << std::setw(2) << std::abs(Seconds());

        if (Microseconds() != 0)
        {
            ss << '.' << std::setw(6) << std::abs(Microseconds());
        }

        return ss.str();
    }

private:
    int64_t mTicks{};

    void CalculateTicks(int days, int hours, int minutes, int seconds, int microseconds)
    {
        mTicks = days * kTICKS_PER_DAY + (hours * 3600LL + minutes * 60LL + seconds) * kTICKS_PER_SECOND +
                 microseconds * kTICKS_PER_MICROSECOND;
    }
};

inline std::ostream& operator<<(std::ostream& strm, const TimeSpan& t)
{
    return strm << t.ToString();
}

inline TimeSpan operator+(const TimeSpan& ts1, const TimeSpan& ts2)
{
    return ts1.Add(ts2);
}

inline TimeSpan operator-(const TimeSpan& ts1, const TimeSpan& ts2)
{
    return ts1.Subtract(ts2);
}

inline bool operator==(const TimeSpan& ts1, const TimeSpan& ts2)
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

} // namespace libsgp4
