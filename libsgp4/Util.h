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


#ifndef UTIL_H_
#define UTIL_H_

#ifdef _WIN32
#define NOMINMAX
#include <Windows.h>
#endif

#include "Globals.h"

#include <iostream>
#include <mutex>
#include <algorithm>
#include <cctype>
#include <ctime>
#include <atomic>
#include <sstream>

namespace Util
{

	template
		<typename T>
		bool FromString(const std::string& str, T& val)
	{
		std::stringstream ss(str);
		return !(ss >> val).fail();
	}

	/*
	 * always positive result
	 * Mod(-3,4)= 1   fmod(-3,4)= -3
	 */
	inline double Mod(const double x, const double y)
	{
		if (y == 0.0)
		{
			return x;
		}

		return x - y * floor(x / y);
	}

	inline double WrapNegPosPI(const double a)
	{
		return Mod(a + kPI, kTWOPI) - kPI;
	}

	inline double WrapTwoPI(const double a)
	{
		return Mod(a, kTWOPI);
	}

	inline double WrapNegPos180(const double a)
	{
		return Mod(a + 180.0, 360.0) - 180.0;
	}

	inline double Wrap360(const double a)
	{
		return Mod(a, 360.0);
	}

	inline double DegreesToRadians(const double degrees)
	{
		return degrees * kPI / 180.0;
	}

	inline double RadiansToDegrees(const double radians)
	{
		return radians * 180.0 / kPI;
	}

	inline double AcTan(const double sinx, const double cosx)
	{
		if (cosx == 0.0)
		{
			if (sinx > 0.0)
			{
				return kPI / 2.0;
			}
			else
			{
				return 3.0 * kPI / 2.0;
			}
		}
		else
		{
			if (cosx > 0.0)
			{
				return atan(sinx / cosx);
			}
			else
			{
				return kPI + atan(sinx / cosx);
			}
		}
	}

	// time since unix epoch
	struct UnixTimestamp
	{
		int64_t seconds;
		int64_t nanoseconds;
	};

    #ifdef _WIN32
	    typedef void (WINAPI *FuncT) (LPFILETIME lpSystemTimeAsFileTime);

	    class WindowsUtcTimeSource
	    {
	    public:
		    inline bool IsPrecise() { return m_isPrecise; };

		    static WindowsUtcTimeSource* Instance();
		    UnixTimestamp GetUtcTime();

	    private:
		    static std::atomic<WindowsUtcTimeSource*> m_instance;
        static std::mutex m_mutex;
		    FuncT m_getUtcTimeFunc;
		    bool m_isPrecise;

		    WindowsUtcTimeSource();
		    WindowsUtcTimeSource(const WindowsUtcTimeSource&) = delete;
		    WindowsUtcTimeSource& operator=(WindowsUtcTimeSource const&) = delete;
	    };
    #endif

	inline UnixTimestamp CurrentTime()
	{
		UnixTimestamp stamp;
    		#ifdef _WIN32
                	stamp = WindowsUtcTimeSource::Instance()->GetUtcTime();		
        	#else
      		  	struct timespec ts;
        		if (clock_gettime(CLOCK_REALTIME, &ts) == 0)
			{
				stamp.seconds = ts.tv_sec;
				stamp.nanoseconds = ts.tv_nsec;
			}
			else
			{
				throw 1;
			}	    
        #endif
		return stamp;
	}

	static inline void TrimLeft(std::string &s)
	{
    		s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch)
		{
	    		return !std::isspace(ch);	
		}));
	}

	static inline void TrimRight(std::string &s) {
		s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch)
		{
			return !std::isspace(ch);
		}).base(), s.end());
	}

	static inline void Trim(std::string &s)
	{
		TrimLeft(s);
		TrimRight(s);
	}

}

#endif
