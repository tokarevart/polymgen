#pragma once


template<class T>
class no_garb_ptr
{
private:
	bool isValid = false;
	T* _ptr;

public:
	//T*& const get() const;

	no_garb_ptr();
	//no_garb_ptr(const T*& ptr); // Need func ptr
	~no_garb_ptr();
};