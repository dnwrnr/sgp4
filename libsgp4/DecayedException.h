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

#include "DateTime.h"
#include "Vector.h"

#include <stdexcept>
#include <string>

namespace libsgp4
{

/**
 * @brief The exception that the SGP4 class throws when a satellite decays.
 */
class DecayedException : public std::runtime_error
{
public:
    /**
     * Constructor
     * @param[in] dt time of the event
     * @param[in] pos position of the satellite at dt
     * @param[in] vel velocity of the satellite at dt
     */
    DecayedException(const DateTime& dt, const Vector& pos, const Vector& vel)
        : runtime_error("Satellite decayed")
        , mDt(dt)
        , mPos(pos)
        , mVel(vel)
    {
    }

    /**
     * @returns the date
     */
    DateTime Decayed() const
    {
        return mDt;
    }

    /**
     * @returns the position
     */
    Vector Position() const
    {
        return mPos;
    }

    /**
     * @returns the velocity
     */
    Vector Velocity() const
    {
        return mVel;
    }

private:
    DateTime mDt;
    Vector mPos;
    Vector mVel;
};

} // namespace libsgp4
