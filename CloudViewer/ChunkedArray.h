#ifndef CSCHUNKEDARRAY_H
#define CSCHUNKEDARRAY_H

#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <vector>
#include <assert.h>

#define CHUNK_INDEX_BIT_DEC 16
static const unsigned int MAX_NUMBER_OF_ELEMENTS_PER_CHUNK = (1 << CHUNK_INDEX_BIT_DEC);// 2^16
static const unsigned int ELEMENT_INDEX_BIT_MASK = MAX_NUMBER_OF_ELEMENTS_PER_CHUNK - 1;

//! template class of generic chunked array
/** split a large array into many smaller arrays(called chunked array)
\ to avoid memory limit of continuous memory allocation
\ D is the dimension of the array(D > 0)
\ Assume that each dimension of the array has the same amount of elements
**/
template <unsigned D, typename ElementType> class ChunkedArray
{
public:
	ChunkedArray()
	{

	}
	virtual ~ChunkedArray();

	bool reserve(unsigned int elementCount)
	{

		return true;
	}

	inline unsigned int size() { return element_count_; }

	inline void addElement(const ElementType* element)
	{

	}

	inline void setValue(unsigned int index, const ElementType* element)
	{
		assert(index < element_count_);
	}

	void clear();


private:
	unsigned int element_count_;
	std::vector<ElementType*> chunks_;
	ElementType min_vals_[D];
	ElementType max_vals_[D];
}


#endif // !CSCHUNKEDARRAY
