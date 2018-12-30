#pragma once
#include <string>
#include <chrono>
#include <vector>
#include <fstream>

enum TimeScale
{
	microseconds,
	milliseconds,
	seconds
};

class Timer
{
	bool _isStarted = false;
	bool _isPaused  = false;

	std::chrono::high_resolution_clock::time_point _startTime;
	std::chrono::high_resolution_clock::time_point _pauseBegin;

	std::chrono::high_resolution_clock::duration              _pauseDuration = std::chrono::high_resolution_clock::duration::zero();
	std::vector<std::chrono::high_resolution_clock::duration> _durations;
	std::vector<std::string>                                  _descriptions;
	std::vector<size_t>                                       _nestingLevels;

public:
	std::string globalDescription;

	// If measurement already started it restarts.
	void Start   (const std::string& description = "", const size_t& nestingLevel = 0);
	void Pause   ();
	void Continue();
	void Stop    ();

	size_t Elapsed                             (TimeScale ts) const;
	size_t GetLastDuration                     (TimeScale ts) const;
	size_t GetDuration    (const size_t& index, TimeScale ts) const;

	void GenerateDurationsLogFile(std::string filename = "time_log.txt") const;

	void   ClearDurationsData ();
	size_t GetDurationDataSize() const;

	Timer();
	Timer(const std::string& description);
	~Timer();
};