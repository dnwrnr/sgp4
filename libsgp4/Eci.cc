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


#include "Eci.h"

#include "Globals.h"
#include "Util.h"

namespace libsgp4
{

/**
 * Converts a DateTime and Geodetic position to Eci coordinates
 * @param[in] dt the date
 * @param[in] geo the geodetic position
 */
void Eci::ToEci(const DateTime& dt, const CoordGeodetic &geo)
{
    /*
     * set date
     */
    mDt = dt;

    static const double mfactor = kTWOPI * (kOMEGA_E / kSECONDS_PER_DAY);

    /*
     * Calculate Local Mean Sidereal Time for observers longitude
     */
    const double theta = mDt.ToLocalMeanSiderealTime(geo.longitude);

    /*
     * take into account earth flattening
     */
    const double c = 1.0
        / sqrt(1.0 + kF * (kF - 2.0) * pow(sin(geo.latitude), 2.0));
    const double s = pow(1.0 - kF, 2.0) * c;
    const double achcp = (kXKMPER * c + geo.altitude) * cos(geo.latitude);

    /*
     * X position in km
     * Y position in km
     * Z position in km
     * W magnitude in km
     */
    mPosition.x = achcp * cos(theta);
    mPosition.y = achcp * sin(theta);
    mPosition.z = (kXKMPER * s + geo.altitude) * sin(geo.latitude);
    mPosition.w = mPosition.Magnitude();

    /*
     * X velocity in km/s
     * Y velocity in km/s
     * Z velocity in km/s
     * W magnitude in km/s
     */
    mVelocity.x = -mfactor * mPosition.y;
    mVelocity.y = mfactor * mPosition.x;
    mVelocity.z = 0.0;
    mVelocity.w = mVelocity.Magnitude();
}

/**
 * @returns the position in geodetic form
 */
CoordGeodetic Eci::ToGeodetic() const
{
    const double theta = Util::AcTan(mPosition.y, mPosition.x);

    const double lon = Util::WrapNegPosPI(theta
            - mDt.ToGreenwichSiderealTime());

    const double r = sqrt((mPosition.x * mPosition.x)
            + (mPosition.y * mPosition.y));

    static const double e2 = kF * (2.0 - kF);

    double lat = Util::AcTan(mPosition.z, r);
    double phi = 0.0;
    double c = 0.0;
    int cnt = 0;

    do
    {
        phi = lat;
        const double sinphi = sin(phi);
        c = 1.0 / sqrt(1.0 - e2 * sinphi * sinphi);
        lat = Util::AcTan(mPosition.z + kXKMPER * c * e2 * sinphi, r);
        cnt++;
    }
    while (fabs(lat - phi) >= 1e-10 && cnt < 10);

    const double alt = r / cos(lat) - kXKMPER * c;

    return CoordGeodetic(lat, lon, alt, true);
}

} // namespace libsgp4
