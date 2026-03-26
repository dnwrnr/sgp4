#include <SGP4.h>
#include <Tle.h>
#include <Vector.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <Globals.h>

using namespace libsgp4;

int main() {
    std::string l1 = "1 18960U 88022A   20001.30746459  .00000074  00000-0  16196-4 0  9998";
    std::string l2 = "2 18960  99.0194  66.6961 0011855 238.9327 121.1118 14.13525141129339";
    Tle tle("18960", l1, l2);
    SGP4 model(tle);
    
    double nodeo = model.elements_.AscendingNode();
    double inclo = model.elements_.Inclination();
    double omgo = model.elements_.ArgumentPerigee();
    double mo = model.elements_.MeanAnomoly();
    double r = model.elements_.RecoveredSemiMajorAxis();

    double snn = sin(nodeo), csn = cos(nodeo);
    double sni = sin(inclo), csi = cos(inclo);

    double u = omgo + mo;
    double snu = sin(u), csu = cos(u);

    // Orientation vectors from SGP4.cc (580+)
    double xmx = -snn * csi;
    double xmy = csn * csi;
    double ux = xmx * snu + csn * csu;
    double uy = xmy * snu + snn * csu;
    double uz = sni * snu;
    
    std::cout << "Original UX: " << ux << " UY: " << uy << " UZ: " << uz << std::endl;
    
    // Position using full rotation Rz(node)*Rx(inc)*Rz(argp) = [r*cos(u), r*sin(u), 0]
    // My previous "manual" rotation Rz(node)*Rx(inc)*Rz(u) from [1, 0, 0] matched this.
    
    // What if inclination is reversed?
    double ux_r = csn * csu - snn * csi * snu;
    double uy_r = snn * csu + csn * csi * snu;
    double uz_r = sni * snu;
    std::cout << "ux_r: " << ux_r << " uy_r: " << uy_r << " uz_r: " << uz_r << std::endl;

    // What if Z is swapped with Y?
    // What if the matrix is Rz(node) * Ry(inc)?
    // Or Rx(inc) * Rz(node)?
    
    Vector p0 = model.FindPosition(0).Position();
    double r0 = sqrt(p0.x*p0.x + p0.y*p0.y + p0.z*p0.z);
    std::cout << "Scalar:      " << p0.x / r0 << " " << p0.y / r0 << " " << p0.z / r0 << std::endl;

    return 0;
}
