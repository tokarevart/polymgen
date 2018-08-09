#pragma once
#include <memory>

template<class T>
class unique_ptr_helper
{
public:
	std::unique_ptr<T> _uniquePtr;

	std::unique_ptr<T>* GetPtrToUniquePtr()
	{
		return &_uniquePtr;
	}

	unique_ptr_helper() {}
	unique_ptr_helper(T* const& ptr)
	{
		_uniquePtr.reset(ptr);
	}
	~unique_ptr_helper()
	{
		_uniquePtr.release();
	}
};