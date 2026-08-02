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

#include "CoordGeodetic.h"
#include "Eci.h"

namespace libsgp4
{

class DateTime;
struct CoordTopocentric;

/**
 * @brief Stores an observers location in Eci coordinates.
 */
class Observer
{
public:
    /**
     * Constructor
     * @param[in] latitude observers latitude in degrees
     * @param[in] longitude observers longitude in degrees
     * @param[in] altitude observers altitude in kilometers
     */
    Observer(double latitude,
            double longitude,
            double altitude)
        : mGeo(latitude, longitude, altitude)
        , mEci(DateTime(), mGeo)
    {
    }

    /**
     * Constructor
     * @param[in] geo the observers position
     */
    explicit Observer(const CoordGeodetic &geo)
        : mGeo(geo)
        , mEci(DateTime(), geo)
    {
    }

    /**
     * Set the observers location
     * @param[in] geo the observers position
     */
    void SetLocation(const CoordGeodetic& geo)
    {
        mGeo = geo;
        mEci.Update(mEci.GetDateTime(), mGeo);
    }

    /**
     * Get the observers location
     * @returns the observers position
     */
    CoordGeodetic GetLocation() const
    {
        return mGeo;
    }

    /**
     * Get the look angle for the observers position to the object
     * @param[in] eci the object to find the look angle to
     * @returns the lookup angle
     */
    CoordTopocentric GetLookAngle(const Eci &eci);

private:
    /**
     * @param[in] dt the date to update the observers position for
     */
    void Update(const DateTime &dt)
    {
        if (mEci != dt)
        {
            mEci.Update(dt, mGeo);
        }
    }

    /** the observers position */
    CoordGeodetic mGeo;
    /** the observers Eci for a particular time */
    Eci mEci;
};

} // namespace libsgp4
