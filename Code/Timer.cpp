#include "Timer.h"

size_t DurationCast(std::chrono::high_resolution_clock::duration duration, TimeScale timeScale)
{
    switch (timeScale)
    {
    case TimeScale::microseconds:
        return std::chrono::duration_cast<std::chrono::microseconds>(duration).count();

    case TimeScale::milliseconds:
        return std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();

    case TimeScale::seconds:
        return std::chrono::duration_cast<std::chrono::seconds>(duration).count();
        
    default:
        throw std::invalid_argument("Wrong time scale.");
    }
}

void Timer::Start(const std::string& description, const size_t& nestingLevel)
{
    if (_isStarted)
        Stop();

    _descriptions.push_back(description);
    _nestingLevels.push_back(nestingLevel);

    _isStarted = true;
    _startTime = std::chrono::high_resolution_clock::now();
    _pauseDuration = std::chrono::high_resolution_clock::duration::zero();

}

void Timer::Pause()
{
    _isPaused = true;
    _pauseBegin = std::chrono::high_resolution_clock::now();
}

void Timer::Continue()
{
    if (!_isPaused)
        return;

    _isPaused = false;
    _pauseDuration += std::chrono::high_resolution_clock::now() - _pauseBegin;
}

size_t Timer::Elapsed(TimeScale timeScale) const
{
    return DurationCast(
          (std::chrono::high_resolution_clock::now() - _startTime) 
        - (_pauseDuration + std::chrono::high_resolution_clock::now() - _pauseBegin), 
        timeScale);
}

void Timer::Stop()
{
    if (!_isStarted)
        return;

    Continue();

    _isStarted = false;
    _durations.push_back((std::chrono::high_resolution_clock::now() - _startTime) - _pauseDuration);
}

size_t Timer::GetLastDuration(TimeScale timeScale) const
{
    return DurationCast(_durations.back(), timeScale);
}

size_t Timer::GetDuration(const size_t& index, TimeScale timeScale) const
{    
    return DurationCast(_durations.at(index), timeScale);
}

void Timer::GenerateDurationsLogFile(std::string filename) const
{
    std::ofstream file(filename, std::ios::app);
    
    size_t max_descriptions_size = 0ull;
    for (size_t i = 0ull, buf, max = _descriptions.size(); i < max; i++)
    {
        buf = _descriptions[i].size();
        if (buf > max_descriptions_size)
        {
            max_descriptions_size = buf;
        }
    }

    size_t extended_line_description_size = max_descriptions_size + 10ull;

    file << globalDescription;
    std::string buf_line;
    for (size_t i = 0, max = _durations.size(); i < max; i++)
    {
        file << '\n';
        buf_line = 
              ">>  " 
            + std::string(_nestingLevels[i] * 4ull, ' ') 
            + _descriptions[i] 
            + ' ' 
            + std::string(extended_line_description_size - 1ull - _descriptions[i].size(), '-') 
            + ">>  ";
        file << buf_line << (double)GetDuration(i, microseconds) * 1e-6 << " s";
    }
    file << "\n\n";
    size_t total_dur_mcs = 0ull;
    for (size_t i = 0ull, max = _durations.size(); i < max; i++)
    {
        total_dur_mcs += GetDuration(i, microseconds);
    }
    file << ">>  Total time: " << (double)total_dur_mcs * 1e-6 << " s";
}

void Timer::ClearDurationsData()
{
    _durations.clear();
    _durations.shrink_to_fit();
}

size_t Timer::GetDurationDataSize() const
{
    return _durations.size();
}

Timer::Timer() {}

Timer::Timer(const std::string& description)
{
    globalDescription = description;
}

Timer::~Timer() {}