#include <libsgp4/CsvTleLoader.h>
#include <libsgp4/SGP4.h>
#include <libsgp4/Eci.h>
#include <filesystem>
#include <iostream>
#include <iomanip>

int main()
{
    std::filesystem::path csv_path =
        std::filesystem::path(CSV_DATA_DIR) / "geodetic.csv";

    std::vector<libsgp4::Tle> tles = libsgp4::LoadCsvTleFile(csv_path.string());
    if (tles.empty())
    {
        std::cerr << "No TLEs loaded from " << csv_path << std::endl;
        return 1;
    }

    const libsgp4::Tle& tle = tles.front();
    libsgp4::SGP4 sgp4(tle);

    libsgp4::DateTime start(2026, 7, 27, 0, 0, 0);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Satellite: " << tle.Name()
              << " (NORAD " << tle.NoradNumber() << ")" << std::endl;
    std::cout << "Start:     " << start << std::endl;
    std::cout << "Steps:     10 x 15 min" << std::endl;
    std::cout << std::string(72, '-') << std::endl;

    for (int i = 0; i <= 10; ++i)
    {
        libsgp4::DateTime dt = start.AddMinutes(i * 15);
        libsgp4::Eci eci = sgp4.FindPosition(dt);
        libsgp4::Vector pos = eci.Position();
        libsgp4::Vector vel = eci.Velocity();

        std::cout << dt
                  << "  pos=[" << pos.x << ", " << pos.y << ", " << pos.z << "]"
                  << "  vel=[" << vel.x << ", " << vel.y << ", " << vel.z << "]"
                  << std::endl;
    }

    return 0;
}
