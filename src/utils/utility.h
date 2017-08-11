#pragma once
#include "../pch.h"

typedef unsigned char uchar;


namespace wcv {
#define MALLOC_ALIGN 16
#define PI   3.1415926535897932384626433832795
#define _2PI	 6.283185307179586476925286766559

	template<typename _Tp>
	/**@brief alligned memory*/
	_Tp* alignPtr(_Tp* ptr, int align = (int)sizeof(_Tp)) {
		assert((align & (align - 1)) == 0);
		return (_Tp*)(((size_t)ptr + align - 1) & -align);
	}

	void* fastAlloc(size_t size);
	void fastFree(void* _mem);
	void * aligned_malloc(size_t size, int align);
	void aligned_free(void * _mm);

};