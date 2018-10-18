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

#include <Util.h>

#ifdef _WIN32

namespace Util
{

	std::atomic<WindowsUtcTimeSource*> WindowsUtcTimeSource::m_instance = std::atomic<WindowsUtcTimeSource*>(nullptr);
	std::mutex WindowsUtcTimeSource::m_mutex;

	WindowsUtcTimeSource* WindowsUtcTimeSource::Instance()
	{
		WindowsUtcTimeSource* tmp = m_instance.load(std::memory_order_relaxed);
		std::atomic_thread_fence(std::memory_order_acquire);
		if (tmp == nullptr)
		{
			std::lock_guard<std::mutex> lock(m_mutex);
			if (tmp == nullptr)
			{
				tmp = new WindowsUtcTimeSource();
				std::atomic_thread_fence(std::memory_order_release);
				m_instance.store(tmp, std::memory_order_relaxed);
			}

		}
		return tmp;
	}

	WindowsUtcTimeSource::WindowsUtcTimeSource()
	{
		// we check to see if GetSystemTimePreciseAsFileTime exists. If it does, use that. Otherwise, use the non-precise version
		HINSTANCE hDLL = LoadLibrary("Kernel32.dll");
		m_getUtcTimeFunc = (FuncT)GetProcAddress((HMODULE)hDLL, "GetSystemTimePreciseAsFileTime");

		if (m_getUtcTimeFunc != nullptr)
		{
			m_isPrecise = true;
		}
		else
		{
			m_isPrecise = false;
			m_getUtcTimeFunc = (FuncT)GetSystemTimeAsFileTime;
		}
	};

	UnixTimestamp WindowsUtcTimeSource::GetUtcTime()
	{
		UnixTimestamp stamp;
		int ticks_per_second = 10000000; // there are 1e7 100ns ticks in a second
		int nanoseconds_per_tick = 100;
		ULONGLONG seconds_from_windows_to_unix_epoch = 11644473600; // seconds from 1601-01-01T00:00:00Z to 1970-01-01T00:00:00Z
		FILETIME ft;
		m_getUtcTimeFunc(&ft);
		ULARGE_INTEGER ui;
		ui.LowPart = ft.dwLowDateTime;
		ui.HighPart = ft.dwHighDateTime;
		int64_t diff = ui.QuadPart - seconds_from_windows_to_unix_epoch * ticks_per_second;
		stamp.seconds = diff / ticks_per_second;
		stamp.nanoseconds = (diff - stamp.seconds * ticks_per_second) * nanoseconds_per_tick;
		return stamp;
	}

}

#endif
