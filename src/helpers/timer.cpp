#include "helpers/timer.h"
#include <stdexcept>

using namespace tva;




size_t Timer::durationCast(hi_res_clock::duration duration, TimeScale timeScale) const
{
    switch (timeScale)
    {
    case TimeScale::Microseconds:
        return static_cast<size_t>(std::chrono::duration_cast<std::chrono::microseconds>(duration).count());

    case TimeScale::Milliseconds:
        return static_cast<size_t>(std::chrono::duration_cast<std::chrono::milliseconds>(duration).count());

    case TimeScale::Seconds:
        return static_cast<size_t>(std::chrono::duration_cast<std::chrono::seconds>(duration).count());
    }
    throw std::exception();
}




void Timer::start(size_t nestingLevel)
{
    if (m_isStarted)
        stop();

    m_nestingLevels.push_back(nestingLevel);

    m_isStarted = true;
    m_startTime = hi_res_clock::now();
    m_pauseDuration = hi_res_clock::duration::zero();

}


void Timer::suspend()
{
    m_isSuspended = true;
    m_pauseBegin = hi_res_clock::now();
}


void Timer::resume()
{
    if (!m_isSuspended)
        return;

    m_isSuspended = false;
    m_pauseDuration += hi_res_clock::now() - m_pauseBegin;
}


void Timer::stop()
{
    if (!m_isStarted)
        return;

    resume();

    m_isStarted = false;
    m_durations.push_back((hi_res_clock::now() - m_startTime) - m_pauseDuration);
}




size_t Timer::elapsed(TimeScale timeScale) const
{
    return durationCast(
        (hi_res_clock::now() - m_startTime)
        - (m_pauseDuration + hi_res_clock::now() - m_pauseBegin),
        timeScale);
}


size_t Timer::getLastDuration(TimeScale timeScale) const
{
    return durationCast(m_durations.back(), timeScale);
}


size_t Timer::getDuration(size_t index, TimeScale timeScale) const
{    
    return durationCast(m_durations.at(index), timeScale);
}




void Timer::clearDurationsData()
{
    m_durations.clear();
    m_durations.shrink_to_fit();
}


size_t Timer::getDurationsNum() const
{
    return m_durations.size();
}




Timer::Timer() {}


Timer::~Timer() {}
