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
template <unsigned int D, typename ElementType> 
class ChunkedArray
{
public:
	ChunkedArray()
		:element_count_(0), capcity_(0)
	{
		memset(min_vals_, 0, sizeof(ElementType)*D);
		memset(min_vals_, 0, sizeof(ElementType)*D);
	}

	virtual ~ChunkedArray()
	{
		while (!chunks_.empty())
		{
			free(chunks_.back());
			chunks_.pop_back();
		}
	}


	inline unsigned int size() 
	{
		return element_count_;
	}

	inline unsigned int capacity()
	{
		return capcity_;
	}

	inline bool addElement(const ElementType* element)
	{
		if (element_count_ < capcity_)
		{
			++element_count_;
			setValue(element_count_-1, element);
			return true;
		}
		else
		{
			ElementType* newChunkMem = (ElementType*)calloc(MAX_NUMBER_OF_ELEMENTS_PER_CHUNK, sizeof(ElementType)*D);
			if (!newChunkMem)
			{
				return false;
			}
			else
			{
				chunks_.push_back(newChunkMem);
				capcity_ += MAX_NUMBER_OF_ELEMENTS_PER_CHUNK;
				++element_count_;
				setValue(element_count_-1, element);
				return true;
			}
		}
	}

	inline void setValue(unsigned int index, const ElementType* v)
	{
		assert(index < element_count_);
		memcpy(value(index), v, D*sizeof(ElementType));

		if (element_count_ == 1)
		{
			memcpy(chunks_[0], min_vals_, D*sizeof(ElementType));
			memcpy(chunks_[0], max_vals_, D*sizeof(ElementType));
		}
		else
		{
			for (int i = 0; i < D; i++)
			{
				if (v[i]<min_vals_[i])
				{
					min_vals_[i] = v[i];
				}else if (v[i]>max_vals_[i])
				{
					max_vals_[i] = v[i];
				}
			}
		}

	}

	inline ElementType* value(unsigned int index)
	{
		assert(index < element_count_);
		return chunks_[index >> CHUNK_INDEX_BIT_DEC] + ((index&ELEMENT_INDEX_BIT_MASK)*D);
	}

	inline ElementType minVal(unsigned int d)
	{
		assert(d < D);
		return min_vals_[d];
	}

	inline ElementType maxVal(unsigned int d)
	{
		assert(d < D);
		return max_vals_[d];
	}

	inline ElementType* minValsPtr()
	{
		return min_vals_;
	}

	inline ElementType* maxValsPtr()
	{
		return max_vals_;
	}


	inline unsigned int dimension()
	{
		return D;
	}

	inline unsigned int chunksCount()
	{
		return (unsigned int)chunks_.size();
	}

	inline ElementType* chunkStartPtr(unsigned int index)
	{
		assert(index < chunks_.size());
		return chunks_[index];
	}

	
	void clear()
	{
		while (!chunks_.empty())
		{
			free(chunks_.back());
			chunks_.pop_back();
		}
		element_count_ = 0;
		capacity() = 0;
		memset(min_vals_, 0, sizeof(ElementType)*D);
		memset(max_vals_, 0, sizeof(ElementType)*D);
	}

	inline unsigned int elementCountOfChunk(unsigned int index)
	{
		assert(index < chunks_.size());
		unsigned int size = (unsigned int)chunks_.size();
		if (index == size-1)
		{
			//the last chunk may not be full
			return (element_count_&ELEMENT_INDEX_BIT_MASK);
		}
		else
		{
			return MAX_NUMBER_OF_ELEMENTS_PER_CHUNK;
		}
	}


private:
	unsigned int element_count_;
	unsigned int capcity_;
	std::vector<ElementType*> chunks_;
	ElementType min_vals_[D];
	ElementType max_vals_[D];
};


#endif // !CSCHUNKEDARRAY
