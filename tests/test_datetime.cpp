#include <gtest/gtest.h>
#include <libsgp4/DateTime.h>

using namespace libsgp4;

TEST(DateTimeConstruction, Default)
{
    DateTime dt;
    EXPECT_EQ(dt.Year(), 1);
    EXPECT_EQ(dt.Month(), 1);
    EXPECT_EQ(dt.Day(), 1);
    EXPECT_EQ(dt.Hour(), 0);
    EXPECT_EQ(dt.Minute(), 0);
    EXPECT_EQ(dt.Second(), 0);
    EXPECT_EQ(dt.Microsecond(), 0);
}

TEST(DateTimeConstruction, YearMonthDay)
{
    DateTime dt(2024, 6, 15);
    EXPECT_EQ(dt.Year(), 2024);
    EXPECT_EQ(dt.Month(), 6);
    EXPECT_EQ(dt.Day(), 15);
}

TEST(DateTimeConstruction, FullDateTime)
{
    DateTime dt(2024, 6, 15, 12, 30, 45, 500000);
    EXPECT_EQ(dt.Year(), 2024);
    EXPECT_EQ(dt.Month(), 6);
    EXPECT_EQ(dt.Day(), 15);
    EXPECT_EQ(dt.Hour(), 12);
    EXPECT_EQ(dt.Minute(), 30);
    EXPECT_EQ(dt.Second(), 45);
    EXPECT_EQ(dt.Microsecond(), 500000);
}

TEST(DateTimeConstruction, YearDayOfYear)
{
    DateTime dt(2024u, 167.0);
    EXPECT_EQ(dt.Year(), 2024);
    EXPECT_EQ(dt.Month(), 6);
    EXPECT_EQ(dt.Day(), 15);
}

TEST(DateTimeIsLeapYear, KnownLeapYears)
{
    EXPECT_TRUE(DateTime::IsLeapYear(2000));
    EXPECT_TRUE(DateTime::IsLeapYear(2024));
    EXPECT_TRUE(DateTime::IsLeapYear(1996));
    EXPECT_TRUE(DateTime::IsLeapYear(4));
}

TEST(DateTimeIsLeapYear, NonLeapYears)
{
    EXPECT_FALSE(DateTime::IsLeapYear(1900));
    EXPECT_FALSE(DateTime::IsLeapYear(2023));
    EXPECT_FALSE(DateTime::IsLeapYear(2100));
    EXPECT_FALSE(DateTime::IsLeapYear(1));
}

TEST(DateTimeIsLeapYear, DivisibleBy400)
{
    EXPECT_TRUE(DateTime::IsLeapYear(400));
    EXPECT_TRUE(DateTime::IsLeapYear(800));
    EXPECT_TRUE(DateTime::IsLeapYear(2000));
}

TEST(DateTimeDaysInMonth, CommonYear)
{
    EXPECT_EQ(DateTime::DaysInMonth(2023, 1), 31);
    EXPECT_EQ(DateTime::DaysInMonth(2023, 2), 28);
    EXPECT_EQ(DateTime::DaysInMonth(2023, 3), 31);
    EXPECT_EQ(DateTime::DaysInMonth(2023, 4), 30);
    EXPECT_EQ(DateTime::DaysInMonth(2023, 5), 31);
    EXPECT_EQ(DateTime::DaysInMonth(2023, 6), 30);
    EXPECT_EQ(DateTime::DaysInMonth(2023, 7), 31);
    EXPECT_EQ(DateTime::DaysInMonth(2023, 8), 31);
    EXPECT_EQ(DateTime::DaysInMonth(2023, 9), 30);
    EXPECT_EQ(DateTime::DaysInMonth(2023, 10), 31);
    EXPECT_EQ(DateTime::DaysInMonth(2023, 11), 30);
    EXPECT_EQ(DateTime::DaysInMonth(2023, 12), 31);
}

TEST(DateTimeDaysInMonth, LeapYear)
{
    EXPECT_EQ(DateTime::DaysInMonth(2024, 2), 29);
    EXPECT_EQ(DateTime::DaysInMonth(2000, 2), 29);
}

TEST(DateTimeIsValidYearMonthDay, ValidDates)
{
    EXPECT_TRUE(DateTime::IsValidYearMonthDay(2024, 6, 15));
    EXPECT_TRUE(DateTime::IsValidYearMonthDay(1, 1, 1));
    EXPECT_TRUE(DateTime::IsValidYearMonthDay(9999, 12, 31));
}

TEST(DateTimeIsValidYearMonthDay, InvalidDates)
{
    EXPECT_FALSE(DateTime::IsValidYearMonthDay(0, 1, 1));
    EXPECT_FALSE(DateTime::IsValidYearMonthDay(10000, 1, 1));
    EXPECT_FALSE(DateTime::IsValidYearMonthDay(2024, 0, 1));
    EXPECT_FALSE(DateTime::IsValidYearMonthDay(2024, 13, 1));
    EXPECT_FALSE(DateTime::IsValidYearMonthDay(2024, 1, 0));
    EXPECT_FALSE(DateTime::IsValidYearMonthDay(2024, 1, 32));
    EXPECT_FALSE(DateTime::IsValidYearMonthDay(2023, 2, 29));
}

TEST(DateTimeDayOfYear, KnownValues)
{
    DateTime dt(2024, 1, 1);
    EXPECT_EQ(dt.DayOfYear(2024, 1, 1), 1);

    DateTime dt2(2024, 12, 31);
    EXPECT_EQ(dt2.DayOfYear(2024, 12, 31), 366);

    DateTime dt3(2023, 12, 31);
    EXPECT_EQ(dt3.DayOfYear(2023, 12, 31), 365);
}

TEST(DateTimeDayOfWeek, KnownValues)
{
    DateTime dt(2024, 6, 15);
    EXPECT_EQ(dt.DayOfWeek(), 6);

    DateTime dt2(2024, 6, 14);
    EXPECT_EQ(dt2.DayOfWeek(), 5);

    DateTime dt3(2024, 6, 16);
    EXPECT_EQ(dt3.DayOfWeek(), 0);
}

TEST(DateTimeArithmetic, AddDays)
{
    DateTime dt(2024, 1, 1);
    DateTime result = dt.AddDays(1.0);
    EXPECT_EQ(result.Year(), 2024);
    EXPECT_EQ(result.Month(), 1);
    EXPECT_EQ(result.Day(), 2);
}

TEST(DateTimeArithmetic, AddSeconds)
{
    DateTime dt(2024, 1, 1, 0, 0, 0);
    DateTime result = dt.AddSeconds(3661.0);
    EXPECT_EQ(result.Hour(), 1);
    EXPECT_EQ(result.Minute(), 1);
    EXPECT_EQ(result.Second(), 1);
}

TEST(DateTimeArithmetic, AddMonths)
{
    DateTime dt(2024, 1, 15);
    DateTime result = dt.AddMonths(1);
    EXPECT_EQ(result.Year(), 2024);
    EXPECT_EQ(result.Month(), 2);
    EXPECT_EQ(result.Day(), 15);
}

TEST(DateTimeArithmetic, AddMonthsEndOfMonth)
{
    DateTime dt(2024, 1, 31);
    DateTime result = dt.AddMonths(1);
    EXPECT_EQ(result.Year(), 2024);
    EXPECT_EQ(result.Month(), 2);
    EXPECT_EQ(result.Day(), 29);
}

TEST(DateTimeArithmetic, AddYears)
{
    DateTime dt(2024, 2, 29);
    DateTime result = dt.AddYears(1);
    EXPECT_EQ(result.Year(), 2025);
    EXPECT_EQ(result.Month(), 2);
    EXPECT_EQ(result.Day(), 28);
}

TEST(DateTimeJulian, KnownValues)
{
    DateTime dt(2000, 1, 1, 12, 0, 0);
    double julian = dt.ToJulian();
    EXPECT_NEAR(julian, 2451545.0, 1e-6);
}

TEST(DateTimeJ1900, KnownValues)
{
    DateTime dt(2000, 1, 1, 12, 0, 0);
    double j1900 = dt.ToJ1900();
    EXPECT_NEAR(j1900, 36525.0, 1e-6);
}

TEST(DateTimeJ2000, KnownValues)
{
    DateTime dt(2000, 1, 1, 12, 0, 0);
    double j2000 = dt.ToJ2000();
    EXPECT_NEAR(j2000, 0.0, 1e-6);
}

TEST(DateTimeOperators, Equal)
{
    DateTime dt1(2024, 6, 15);
    DateTime dt2(2024, 6, 15);
    EXPECT_TRUE(dt1 == dt2);
}

TEST(DateTimeOperators, NotEqual)
{
    DateTime dt1(2024, 6, 15);
    DateTime dt2(2024, 6, 16);
    EXPECT_TRUE(dt1 != dt2);
}

TEST(DateTimeOperators, LessThan)
{
    DateTime dt1(2024, 6, 15);
    DateTime dt2(2024, 6, 16);
    EXPECT_TRUE(dt1 < dt2);
}

TEST(DateTimeOperators, GreaterThan)
{
    DateTime dt1(2024, 6, 16);
    DateTime dt2(2024, 6, 15);
    EXPECT_TRUE(dt1 > dt2);
}

TEST(DateTimeOperators, AddTimeSpan)
{
    DateTime dt(2024, 1, 1);
    TimeSpan ts(1, 0, 0, 0);
    DateTime result = dt + ts;
    EXPECT_EQ(result.Day(), 2);
}

TEST(DateTimeOperators, SubtractTimeSpan)
{
    DateTime dt(2024, 1, 2);
    TimeSpan ts(1, 0, 0, 0);
    DateTime result = dt - ts;
    EXPECT_EQ(result.Day(), 1);
}

TEST(DateTimeOperators, SubtractDateTimes)
{
    DateTime dt1(2024, 1, 2);
    DateTime dt2(2024, 1, 1);
    TimeSpan result = dt1 - dt2;
    EXPECT_EQ(result.Days(), 1);
}

TEST(DateTimeTimeOfDay, KnownValues)
{
    DateTime dt(2024, 1, 1, 12, 30, 45);
    TimeSpan tod = dt.TimeOfDay();
    EXPECT_EQ(tod.Hours(), 12);
    EXPECT_EQ(tod.Minutes(), 30);
    EXPECT_EQ(tod.Seconds(), 45);
}

TEST(DateTimeToString, Format)
{
    DateTime dt(2024, 6, 15, 12, 30, 45, 123456);
    std::string str = dt.ToString();
    EXPECT_NE(str.find("2024"), std::string::npos);
    EXPECT_NE(str.find("06"), std::string::npos);
    EXPECT_NE(str.find("15"), std::string::npos);
    EXPECT_NE(str.find("12:30:45"), std::string::npos);
}

TEST(DateTimeRoundTrip, FullRange)
{
    for (int month = 1; month <= 12; month++)
    {
        for (int day = 1; day <= DateTime::DaysInMonth(2024, month); day++)
        {
            DateTime dt(2024, month, day, 23, 59, 59, 999999);
            EXPECT_EQ(dt.Year(), 2024);
            EXPECT_EQ(dt.Month(), month);
            EXPECT_EQ(dt.Day(), day);
            EXPECT_EQ(dt.Hour(), 23);
            EXPECT_EQ(dt.Minute(), 59);
            EXPECT_EQ(dt.Second(), 59);
            EXPECT_EQ(dt.Microsecond(), 999999);
        }
    }
}

TEST(DateTimeFromTicks, LeapYearBoundary)
{
    DateTime dt(2000, 2, 29);
    EXPECT_EQ(dt.Year(), 2000);
    EXPECT_EQ(dt.Month(), 2);
    EXPECT_EQ(dt.Day(), 29);
}
