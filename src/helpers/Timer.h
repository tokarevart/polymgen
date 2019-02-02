#pragma once
#include <stddef.h>
#include <chrono>
#include <vector>

using std::chrono::high_resolution_clock;
typedef high_resolution_clock hi_res_clock;


namespace tva
{
    class Timer
    {
    public:
        enum TimeScale
        {
            MICROSECONDS,
            MILLISECONDS,
            SECONDS
        };

        // If the measurement already started it restarts.
        void start(size_t nestingLevel = 0);
        void suspend();
        void resume();
        void stop();

        size_t elapsed(TimeScale ts) const;
        size_t getLastDuration(TimeScale ts) const;
        size_t getDuration(size_t index, TimeScale ts) const;

        void   clearDurationsData();
        size_t getDurationsNum() const;

        Timer();
        ~Timer();


    private:
        bool m_isStarted = false;
        bool m_isSuspended = false;

        hi_res_clock::time_point m_startTime;
        hi_res_clock::time_point m_pauseBegin;

        hi_res_clock::duration              m_pauseDuration = hi_res_clock::duration::zero();
        std::vector<hi_res_clock::duration> m_durations;
        std::vector<size_t>                 m_nestingLevels;

        size_t durationCast(hi_res_clock::duration duration, TimeScale timeScale) const;
    };
}
