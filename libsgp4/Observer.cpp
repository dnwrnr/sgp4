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

#include "Observer.h"

#include "CoordTopocentric.h"
#include "SatelliteException.h"

#include <cmath>

namespace libsgp4
{

/*
 * calculate lookangle between the observer and the passed in Eci object
 */
CoordTopocentric Observer::GetLookAngle(const Eci& eci)
{
    /*
     * update the observers Eci to match the time of the Eci passed in
     * if necessary
     */
    Update(eci.GetDateTime());

    /*
     * calculate differences
     */
    Vector rangeRate = eci.Velocity() - mEci.Velocity();
    Vector range = eci.Position() - mEci.Position();

    if (!std::isfinite(range.x) || !std::isfinite(range.y) || !std::isfinite(range.z))
    {
        throw SatelliteException("Error: (range not finite)");
    }

    range.w = range.Magnitude();

    if (!(range.w > 0.0))
    {
        throw SatelliteException("Error: (range zero or not finite)");
    }

    /*
     * Calculate Local Mean Sidereal Time for observers longitude
     */
    double theta = eci.GetDateTime().ToLocalMeanSiderealTime(mGeo.longitude);

    double sinLat = sin(mGeo.latitude);
    double cosLat = cos(mGeo.latitude);
    double sinTheta = sin(theta);
    double cosTheta = cos(theta);

    double topS = sinLat * cosTheta * range.x + sinLat * sinTheta * range.y - cosLat * range.z;
    double topE = -sinTheta * range.x + cosTheta * range.y;
    double topZ = cosLat * cosTheta * range.x + cosLat * sinTheta * range.y + sinLat * range.z;
    double az = atan(-topE / topS);

    if (topS > 0.0)
    {
        az += kPI;
    }

    if (az < 0.0)
    {
        az += 2.0 * kPI;
    }

    double el = asin(Util::Clamp(topZ / range.w, -1.0, 1.0));
    double rate = range.Dot(rangeRate) / range.w;

    /*
     * azimuth in radians
     * elevation in radians
     * range in km
     * range rate in km/s
     */
    return CoordTopocentric(az, el, range.w, rate);
}

} // namespace libsgp4
