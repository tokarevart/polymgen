#include "Timer.h"




size_t tva::Timer::durationCast(std::chrono::high_resolution_clock::duration duration, TimeScale timeScale) const
{
    switch (timeScale)
    {
    case MICROSECONDS:
        return std::chrono::duration_cast<std::chrono::microseconds>(duration).count();

    case MILLISECONDS:
        return std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();

    case SECONDS:
        return std::chrono::duration_cast<std::chrono::seconds>(duration).count();
        
    default:
        throw std::invalid_argument("Wrong time scale.");
    }
}




void tva::Timer::start(size_t nestingLevel)
{
    if (m_isStarted)
        stop();

    m_nestingLevels.push_back(nestingLevel);

    m_isStarted = true;
    m_startTime = hi_res_clock::now();
    m_pauseDuration = hi_res_clock::duration::zero();

}


void tva::Timer::suspend()
{
    m_isSuspended = true;
    m_pauseBegin = hi_res_clock::now();
}


void tva::Timer::resume()
{
    if (!m_isSuspended)
        return;

    m_isSuspended = false;
    m_pauseDuration += hi_res_clock::now() - m_pauseBegin;
}


void tva::Timer::stop()
{
    if (!m_isStarted)
        return;

    resume();

    m_isStarted = false;
    m_durations.push_back((hi_res_clock::now() - m_startTime) - m_pauseDuration);
}




size_t tva::Timer::elapsed(TimeScale timeScale) const
{
    return durationCast(
        (hi_res_clock::now() - m_startTime)
        - (m_pauseDuration + hi_res_clock::now() - m_pauseBegin),
        timeScale);
}


size_t tva::Timer::getLastDuration(TimeScale timeScale) const
{
    return durationCast(m_durations.back(), timeScale);
}


size_t tva::Timer::getDuration(size_t index, TimeScale timeScale) const
{    
    return durationCast(m_durations.at(index), timeScale);
}




void tva::Timer::clearDurationsData()
{
    m_durations.clear();
    m_durations.shrink_to_fit();
}


size_t tva::Timer::getDurationsNum() const
{
    return m_durations.size();
}




tva::Timer::Timer() {}


tva::Timer::~Timer() {}