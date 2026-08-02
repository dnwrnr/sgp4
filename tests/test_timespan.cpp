#include <gtest/gtest.h>
#include <libsgp4/TimeSpan.h>

using namespace libsgp4;

TEST(TimeSpanConstruction, DefaultTicks)
{
    TimeSpan ts(0LL);
    EXPECT_EQ(ts.Days(), 0);
    EXPECT_EQ(ts.Hours(), 0);
    EXPECT_EQ(ts.Minutes(), 0);
    EXPECT_EQ(ts.Seconds(), 0);
    EXPECT_EQ(ts.Microseconds(), 0);
}

TEST(TimeSpanConstruction, HoursMinutesSeconds)
{
    TimeSpan ts(1, 30, 45);
    EXPECT_EQ(ts.Hours(), 1);
    EXPECT_EQ(ts.Minutes(), 30);
    EXPECT_EQ(ts.Seconds(), 45);
    EXPECT_EQ(ts.TotalSeconds(), 5445.0);
}

TEST(TimeSpanConstruction, DaysHoursMinutesSeconds)
{
    TimeSpan ts(2, 3, 30, 45);
    EXPECT_EQ(ts.Days(), 2);
    EXPECT_EQ(ts.Hours(), 3);
    EXPECT_EQ(ts.Minutes(), 30);
    EXPECT_EQ(ts.Seconds(), 45);
    EXPECT_EQ(ts.TotalDays(), 2.0 + 3.0 / 24.0 + 30.0 / 1440.0 + 45.0 / 86400.0);
}

TEST(TimeSpanConstruction, WithMicroseconds)
{
    TimeSpan ts(0, 0, 0, 1, 500000);
    EXPECT_EQ(ts.Seconds(), 1);
    EXPECT_EQ(ts.Microseconds(), 500000);
}

TEST(TimeSpanComponents, MixedValues)
{
    TimeSpan ts(1, 2, 3, 4, 5);
    EXPECT_EQ(ts.Days(), 1);
    EXPECT_EQ(ts.Hours(), 2);
    EXPECT_EQ(ts.Minutes(), 3);
    EXPECT_EQ(ts.Seconds(), 4);
    EXPECT_EQ(ts.Microseconds(), 5);
}

TEST(TimeSpanComponents, Overflow)
{
    TimeSpan ts(0, 25, 0, 0);
    EXPECT_EQ(ts.Days(), 1);
    EXPECT_EQ(ts.Hours(), 1);
}

TEST(TimeSpanTotal, Days)
{
    TimeSpan ts(1, 0, 0, 0);
    EXPECT_DOUBLE_EQ(ts.TotalDays(), 1.0);
}

TEST(TimeSpanTotal, Hours)
{
    TimeSpan ts(0, 1, 0, 0);
    EXPECT_DOUBLE_EQ(ts.TotalHours(), 1.0);
}

TEST(TimeSpanTotal, Minutes)
{
    TimeSpan ts(0, 0, 1, 0);
    EXPECT_DOUBLE_EQ(ts.TotalMinutes(), 1.0);
}

TEST(TimeSpanTotal, Seconds)
{
    TimeSpan ts(0, 0, 0, 1);
    EXPECT_DOUBLE_EQ(ts.TotalSeconds(), 1.0);
}

TEST(TimeSpanAdd, Positive)
{
    TimeSpan ts1(1, 0, 0, 0);
    TimeSpan ts2(0, 12, 0, 0);
    TimeSpan result = ts1.Add(ts2);
    EXPECT_EQ(result.Days(), 1);
    EXPECT_EQ(result.Hours(), 12);
}

TEST(TimeSpanSubtract, Positive)
{
    TimeSpan ts1(2, 0, 0, 0);
    TimeSpan ts2(1, 0, 0, 0);
    TimeSpan result = ts1.Subtract(ts2);
    EXPECT_EQ(result.Days(), 1);
}

TEST(TimeSpanOperators, Plus)
{
    TimeSpan ts1(1, 0, 0, 0);
    TimeSpan ts2(0, 12, 0, 0);
    TimeSpan result = ts1 + ts2;
    EXPECT_EQ(result.Days(), 1);
    EXPECT_EQ(result.Hours(), 12);
}

TEST(TimeSpanOperators, Minus)
{
    TimeSpan ts1(2, 0, 0, 0);
    TimeSpan ts2(1, 0, 0, 0);
    TimeSpan result = ts1 - ts2;
    EXPECT_EQ(result.Days(), 1);
}

TEST(TimeSpanCompare, Equal)
{
    TimeSpan ts1(1, 0, 0, 0);
    TimeSpan ts2(1, 0, 0, 0);
    EXPECT_TRUE(ts1 == ts2);
    EXPECT_FALSE(ts1 != ts2);
    EXPECT_FALSE(ts1 < ts2);
    EXPECT_FALSE(ts1 > ts2);
    EXPECT_TRUE(ts1 <= ts2);
    EXPECT_TRUE(ts1 >= ts2);
}

TEST(TimeSpanCompare, LessThan)
{
    TimeSpan ts1(1, 0, 0, 0);
    TimeSpan ts2(2, 0, 0, 0);
    EXPECT_TRUE(ts1 < ts2);
    EXPECT_TRUE(ts1 <= ts2);
    EXPECT_FALSE(ts1 > ts2);
    EXPECT_FALSE(ts1 >= ts2);
    EXPECT_TRUE(ts1 != ts2);
}

TEST(TimeSpanCompare, GreaterThan)
{
    TimeSpan ts1(2, 0, 0, 0);
    TimeSpan ts2(1, 0, 0, 0);
    EXPECT_TRUE(ts1 > ts2);
    EXPECT_TRUE(ts1 >= ts2);
    EXPECT_FALSE(ts1 < ts2);
    EXPECT_FALSE(ts1 <= ts2);
    EXPECT_TRUE(ts1 != ts2);
}

TEST(TimeSpanToString, Days)
{
    TimeSpan ts(1, 0, 0, 0);
    EXPECT_EQ(ts.ToString(), "01.00:00:00");
}

TEST(TimeSpanToString, HoursMinutesSeconds)
{
    TimeSpan ts(0, 12, 30, 45);
    EXPECT_EQ(ts.ToString(), "12:30:45");
}

TEST(TimeSpanToString, Negative)
{
    TimeSpan ts(-3600000000LL);
    EXPECT_EQ(ts.ToString(), "-01:00:00");
}

TEST(TimeSpanToString, WithMicroseconds)
{
    TimeSpan ts(0, 0, 0, 1, 500000);
    EXPECT_EQ(ts.ToString(), "00:00:01.500000");
}

TEST(TimeSpanTicks, RoundTrip)
{
    int64_t ticks = 123456789LL;
    TimeSpan ts(ticks);
    EXPECT_EQ(ts.Ticks(), ticks);
}

TEST(TimeSpanNegative, Days)
{
    TimeSpan ts(-1, 0, 0, 0);
    EXPECT_EQ(ts.Days(), -1);
}

TEST(TimeSpanNegative, TotalSeconds)
{
    TimeSpan ts(0, 0, -1, 0);
    EXPECT_DOUBLE_EQ(ts.TotalSeconds(), -60.0);
}
