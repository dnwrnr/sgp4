#include <gtest/gtest.h>
#include <libsgp4/Tle.h>
#include <libsgp4/CsvTleLoader.h>
#include <libsgp4/SGP4.h>
#include <libsgp4/Eci.h>
#include <filesystem>
#include <fstream>
#include <random>

using namespace libsgp4;

namespace
{
class TempFile
{
public:
    explicit TempFile(const std::string& content,
                      const std::string& suffix = ".csv")
    {
        auto dir = std::filesystem::temp_directory_path();
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<uint64_t> dis;
        path_ = dir / ("sgp4_test_" + std::to_string(dis(gen)) + suffix);
        std::ofstream f(path_);
        f << content;
    }

    ~TempFile()
    {
        std::error_code ec;
        std::filesystem::remove(path_, ec);
    }

    TempFile(const TempFile&) = delete;
    TempFile& operator=(const TempFile&) = delete;

    const std::filesystem::path& path() const { return path_; }

private:
    std::filesystem::path path_;
};
}

static const char* CSV_LINE_GPS =
    "GPS BIIR-5  (PRN 22),2000-040A,2026-07-26T04:36:54.986112,2.00558057,"
    ".01202129,54.8474,213.6129,302.8892,63.8849,0,U,26407,999,19073,0,"
    ".81E-6,0";

static const char* CSV_LINE_ISS =
    "ISS (ZARYA),1998-067A,2026-07-26T09:51:36.701280,15.49181604,"
    ".00069412,51.6317,105.9876,339.4731,20.5977,0,U,25544,999,57783,"
    ".20165918E-3,.10751E-3,0";

static const char* CSV_LINE_AMAT_BSTAR =
    "OSCAR 7 (AO-7),1974-089B,2026-07-26T15:43:46.085664,12.53698736,"
    ".00125885,101.9910,221.3137,57.0910,60.8894,0,U,7530,999,36546,"
    ".13291256E-3,-.25E-6,0";

TEST(CsvTleFromCsv, BasicFields)
{
    Tle tle = Tle::FromCsv(CSV_LINE_GPS);
    EXPECT_EQ(tle.Name(), "GPS BIIR-5  (PRN 22)");
    EXPECT_EQ(tle.NoradNumber(), 26407u);
    EXPECT_EQ(tle.IntDesignator(), "2000-040A");
    EXPECT_NEAR(tle.MeanMotion(), 2.00558057, 1e-7);
    EXPECT_NEAR(tle.Eccentricity(), 0.01202129, 1e-7);
    EXPECT_NEAR(tle.Inclination(true), 54.8474, 1e-4);
    EXPECT_NEAR(tle.RightAscendingNode(true), 213.6129, 1e-4);
    EXPECT_NEAR(tle.ArgumentPerigee(true), 302.8892, 1e-4);
    EXPECT_NEAR(tle.MeanAnomaly(true), 63.8849, 1e-4);
    EXPECT_NEAR(tle.BStar(), 0.0, 1e-10);
    EXPECT_NEAR(tle.MeanMotionDt2(), 8.1e-7, 1e-12);
    EXPECT_NEAR(tle.MeanMotionDdt6(), 0.0, 1e-10);
    EXPECT_EQ(tle.OrbitNumber(), 19073u);
}

TEST(CsvTleFromCsv, Epoch)
{
    Tle tle = Tle::FromCsv(CSV_LINE_GPS);
    DateTime epoch = tle.Epoch();
    EXPECT_EQ(epoch.Year(), 2026);
    EXPECT_EQ(epoch.Month(), 7);
    EXPECT_EQ(epoch.Day(), 26);
    EXPECT_EQ(epoch.Hour(), 4);
    EXPECT_EQ(epoch.Minute(), 36);
    EXPECT_EQ(epoch.Second(), 54);
}

TEST(CsvTleFromCsv, IssFields)
{
    Tle tle = Tle::FromCsv(CSV_LINE_ISS);
    EXPECT_EQ(tle.Name(), "ISS (ZARYA)");
    EXPECT_EQ(tle.NoradNumber(), 25544u);
    EXPECT_NEAR(tle.Inclination(true), 51.6317, 1e-4);
    EXPECT_NEAR(tle.BStar(), 0.20165918e-3, 1e-9);
    EXPECT_NEAR(tle.MeanMotionDt2(), 0.10751e-3, 1e-9);
    EXPECT_EQ(tle.OrbitNumber(), 57783u);
}

TEST(CsvTleFromCsv, NegativeBStar)
{
    Tle tle = Tle::FromCsv(CSV_LINE_AMAT_BSTAR);
    EXPECT_NEAR(tle.BStar(), 0.13291256e-3, 1e-10);
    EXPECT_NEAR(tle.MeanMotionDt2(), -0.25e-6, 1e-12);
}

TEST(CsvTleFromCsv, EpochWithoutFractionalSeconds)
{
    Tle tle = Tle::FromCsv(
        "TEST SAT,2000-040A,2026-07-26T04:36:54,2.00558057,"
        ".01202129,54.8474,213.6129,302.8892,63.8849,0,U,26407,999,19073,0,"
        ".81E-6,0");
    DateTime epoch = tle.Epoch();
    EXPECT_EQ(epoch.Year(), 2026);
    EXPECT_EQ(epoch.Month(), 7);
    EXPECT_EQ(epoch.Day(), 26);
    EXPECT_EQ(epoch.Hour(), 4);
    EXPECT_EQ(epoch.Minute(), 36);
    EXPECT_EQ(epoch.Second(), 54);
    EXPECT_EQ(epoch.Microsecond(), 0);
}

TEST(CsvTleFromCsv, Line1Line2Empty)
{
    Tle tle = Tle::FromCsv(CSV_LINE_GPS);
    EXPECT_TRUE(tle.Line1().empty());
    EXPECT_TRUE(tle.Line2().empty());
}

TEST(CsvTleFromCsv, InvalidFieldCount)
{
    EXPECT_THROW(Tle::FromCsv("a,b,c"), TleException);
}

TEST(CsvTleFromCsv, InvalidEpoch)
{
    EXPECT_THROW(
        Tle::FromCsv(
            "NAME,2000-040A,not-a-date,2.0,.01,54.0,213.0,302.0,63.0,"
            "0,U,26407,999,19073,0,.81E-6,0"),
        TleException);
}

TEST(CsvTleLoader, LoadFile)
{
    std::string csv_content =
        "OBJECT_NAME,OBJECT_ID,EPOCH,MEAN_MOTION,ECCENTRICITY,"
        "INCLINATION,RA_OF_ASC_NODE,ARG_OF_PERICENTER,MEAN_ANOMALY,"
        "EPHEMERIS_TYPE,CLASSIFICATION_TYPE,NORAD_CAT_ID,ELEMENT_SET_NO,"
        "REV_AT_EPOCH,BSTAR,MEAN_MOTION_DOT,MEAN_MOTION_DDOT\n"
        "TEST SAT 1,99999A,2026-01-01T00:00:00.000000,15.0,"
        ".001,51.6,100.0,200.0,50.0,0,U,99999,999,100,"
        ".1E-5,.2E-6,0\n"
        "TEST SAT 2,99998A,2026-01-01T12:00:00.000000,14.0,"
        ".002,52.0,101.0,201.0,51.0,0,U,99998,999,200,"
        ".3E-5,.4E-6,0\n";

    TempFile tmp(csv_content);

    std::vector<Tle> tles = LoadCsvTleFile(tmp.path().string());
    ASSERT_EQ(tles.size(), 2u);
    EXPECT_EQ(tles[0].NoradNumber(), 99999u);
    EXPECT_EQ(tles[1].NoradNumber(), 99998u);
    EXPECT_EQ(tles[0].Name(), "TEST SAT 1");
    EXPECT_EQ(tles[1].Name(), "TEST SAT 2");
}

TEST(CsvTleLoader, FileNotFound)
{
    std::vector<Tle> tles = LoadCsvTleFile("/nonexistent/file.csv");
    EXPECT_TRUE(tles.empty());
}

TEST(CsvTleLoader, EmptyFile)
{
    TempFile tmp("");

    std::vector<Tle> tles = LoadCsvTleFile(tmp.path().string());
    EXPECT_TRUE(tles.empty());
}

TEST(CsvTleLoader, HeaderOnly)
{
    TempFile tmp("OBJECT_NAME,OBJECT_ID,EPOCH,MEAN_MOTION\n");

    std::vector<Tle> tles = LoadCsvTleFile(tmp.path().string());
    EXPECT_TRUE(tles.empty());
}

TEST(CsvTleEci, PropagateFromCsv)
{
    Tle tle = Tle::FromCsv(CSV_LINE_ISS);
    SGP4 model(tle);
    Eci eci = model.FindPosition(60.0);
    EXPECT_NE(eci.Position().x, 0.0);
    EXPECT_NE(eci.Position().y, 0.0);
    EXPECT_GT(eci.Position().Magnitude(), 0.0);
}
