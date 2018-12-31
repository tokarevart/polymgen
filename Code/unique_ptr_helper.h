#pragma once
#include <memory>

template<class T>
class unique_ptr_helper
{
	std::unique_ptr<T>* _uniquePtr;

public:
	std::unique_ptr<T>* getPtrToUPtr()
	{
		return _uniquePtr;
	}

	unique_ptr_helper() {}
	unique_ptr_helper(T* ptr)
	{
		_uniquePtr = new std::unique_ptr<T>();
		_uniquePtr->reset(ptr);
	}
	~unique_ptr_helper() {}
};