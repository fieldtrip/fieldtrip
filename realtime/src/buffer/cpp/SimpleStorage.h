/*
 * Copyright (C) 2010, Stefan Klanke
 * Donders Institute for Donders Institute for Brain, Cognition and Behaviour,
 * Centre for Cognitive Neuroimaging, Radboud University Nijmegen,
 * Kapittelweg 29, 6525 EN Nijmegen, The Netherlands
 */
#ifndef __SimpleStorage_h
#define __SimpleStorage_h

/** Simple class for managing a contiguous piece of memory with convenient resizing. 
	Pretty close to what std::vector<char> provides, but this implementation better
	fits a truly typeless void *buffer.
*/
class SimpleStorage {
	public:
	
	/** Constructor. Does not allocate anything, sets capacity and size to 0 */
	SimpleStorage() {
		m_data = NULL;
		m_capacity = m_size = 0;
	}
	
	/** Destructor. Just calls clear() */
	~SimpleStorage() {
		clear();
	}
	
	/** Resizes the buffer. If no memory has been allocated so far, exactly 'size' bytes will
		be allocated. Otherwise, if 'size' is bigger than the current capacity, memory will
		be re-allocated. If 'size' is smaller than or equal to capacity, the memory will
		be kept as it is.
		@param size		Desired buffer size
		@return true	On success (ie., the buffer now holds at least 'size' bytes)
				false 	On memory allocation errors (capacity and size are unchanged)
	*/
	bool resize(unsigned int size) {
		if (size <= m_capacity) {
			m_size = size;
			return true;
		}
		if (m_data != NULL) {
			void *nd = realloc(m_data, size);
			if (nd == NULL) return false;
			m_data = nd;
		} else {
			m_data = malloc(size);
			if (m_data == NULL) return false;
		}
		m_size = m_capacity = size;
		return true;
	}
	
	/** Free any internally allocated memory and set capacity and size to 0 */
	void clear() {
		if (m_data!=NULL) free(m_data);
		m_data = NULL;
		m_size = m_capacity = 0;
	}	
	
	/** Returns pointer to internally allocated memory */
	void *data() { return m_data; }
	
	/** Returns (minimum guaranteed) size of the buffer */
	unsigned int size() { return m_size; }
	
	protected:
	
	unsigned int m_size;
	unsigned int m_capacity;
	void *m_data;
};

#endif
