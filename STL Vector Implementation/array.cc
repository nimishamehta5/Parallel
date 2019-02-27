//
// Created by bpswe on 9/24/2018.
//

//implement your array code here
#include "array.h"
#include <iostream>
using namespace std;

template <typename T>
array<T>::array()
{
	m_elements=nullptr;
	m_size=0;
	m_reserved_size=0;
}

template <typename T>
array<T>::array(std::initializer_list<T> list) : m_size(list.size())
{
	int i;
	m_reserved_size=m_size;
	m_elements=(T*) malloc(sizeof(T)*m_size);
	auto it=list.begin();
	for(i=0;i<m_size;i++)
	{
		new(&m_elements[i]) T(*it);
		it++;
	}
}

template <typename T>
array<T>::array(const array& arr)
{
	int i;
	m_size=arr.m_size;
	m_reserved_size=arr.m_reserved_size;
	m_elements=(T*) malloc(m_reserved_size*sizeof(T));	/*TODO*/
	for(i=0;i<m_size;i++)
	{
		/*new (address) constructor*/
		new (&m_elements[i]) T(arr.m_elements[i]);
		//memcpy
	}
}

template <typename T>
array<T>::array(array&& arr)
{	/*In a move constructor, you are moving data from one block in memory 
	to another bock in memory, so you need to de-allocate/free the original block*/
	m_size=arr.m_size;
	m_reserved_size=arr.m_reserved_size;
	m_elements=arr.m_elements;

	arr.m_elements=nullptr;
	arr.m_size=0;
	arr.m_reserved_size=0;
}

template <typename T>
array<T>::array(size_t sz)
{
	m_size=0;
	m_reserved_size=sz;
	m_elements=(T*)malloc(sz*sizeof(T));
}

template <typename T>
array<T>::array(size_t n, const T& t)
{	
	int i;
	/*m_reserved_size=n;?*/
	m_size=n;
	m_elements=(T*)malloc(n*sizeof(T));

	for(i=0;i<m_size;i++)
	{
		new(&m_elements[i]) T(t);
	}
}

template <typename T>
array<T>::~array()
{
	if(m_elements!=nullptr)
	{
		for(int i=0;i<m_size;i++)
		{
			m_elements[i].~T();
		}
		free(m_elements);
	}
	m_size=0;
	m_reserved_size=0;
//don't delete?
	m_elements=nullptr;
}

template <typename T>
void array<T>::reserve(size_t n)
{	
	int i;
	//if(n>m_reserved_size)
	//{
		T* temp=(T*) malloc(sizeof(T)*n);
		for(i=0;i<m_size;i++)
		{
			new(&temp[i]) T(m_elements[i]);
			m_elements[i].~T();
		}
		free(m_elements);
		m_elements=temp;
		m_reserved_size=n;
	//}
}

template <typename T>
void array<T>::push_back(const T& x)
{	
	int i;
	if(m_size<m_reserved_size)
	{
		new(&m_elements[m_size]) T(x);
		m_size++;
	}
	else
	{
		//creating new array
		m_reserved_size++;
		T* newarr=(T*) malloc(m_reserved_size * sizeof(T));
		for(i=0;i<m_size;i++)
		{
			new(&newarr[i]) T(std::move(m_elements[i])); //std::move
			m_elements[i].~T();
		}
		free(m_elements);
		new(&newarr[m_size]) T(x);
		m_size++;
		m_elements=newarr;
	}

}

template <typename T>
void array<T>::push_front(const T& x)
{	
	int i;
	if(m_size<m_reserved_size)
	{	
		/*downward shift*/
		T* newarr=(T*) malloc(m_reserved_size * sizeof(T));
		for(i=0;i<m_size;i++)
		{
			new(&newarr[i+1]) T(std::move(m_elements[i]));
			m_elements[i].~T();
		}
		free(m_elements);
		new(&newarr[0]) T(x);
		m_elements=newarr;
	}
	else
	{
		/*creating new array and downward shift*/
		m_reserved_size++;
		T* newarr=(T*) malloc(m_reserved_size * sizeof(T));
		for(i=0;i<m_size;i++)
		{
			new(&newarr[i+1]) T(std::move(m_elements[i]));
			m_elements[i].~T();
		}
		free(m_elements);
		new(&newarr[0]) T(x);
		m_elements=newarr;
	}
	m_size++;
}

template <typename T>
void array<T>::pop_back()
{
	if(m_size>0)
	{
		m_elements[m_size-1].~T();
		m_size--;
	}
}

template <typename T>
void array<T>::pop_front()
{
	if(m_size>0)
	{
		for(int i=0;i<m_size-1;i++)
		{
			m_elements[i]=std::move(m_elements[i+1]); //no new. m_el[i]=std::move[i+1]
		}
		pop_back();
	}
}

template <typename T>
T& array<T>::front() const
{
	return m_elements[0];
}

template <typename T>
T& array<T>::back() const
{
	return m_elements[m_size-1];
}

template <typename T>
const T& array<T>::operator[](size_t x) const
{
	return m_elements[x];
}

template <typename T>
T& array<T>::operator[](size_t x)
{
	return m_elements[x];
}

template <typename T>
size_t array<T>::length() const
{
	return m_size;
}

template <typename T>
bool array<T>::empty() const
{
	if(m_size==0)
		{return true;}
	else
		{return false;}
}

template <typename T>
void array<T>::clear()
{	
	int i;
	for(i=0;i<m_size;i++)
	{
		m_elements[i].~T();
	}
	m_size=0;
}

template <typename T>
array_iterator<T> array<T>::begin() const
{
	return(array_iterator<T>(m_elements));
}

template <typename T>
array_iterator<T> array<T>::end() const
{
	return(array_iterator<T>(m_elements+m_size));
}

template <typename T>
void array<T>::erase(const array_iterator<T>& it)
{
	/*use another iterator to find element indicated by 'it', delete it and shift others*/
	int i;
	int pos=0;
	array_iterator<T> ti;
	ti=begin();
	while(ti!=end())
	{
		if(ti == it)
			{break;}
		pos++;
		ti++;
	}
	/*destroy m_el[pos]*/
	m_elements[pos].~T();
	for(i=pos;i<m_size;i++)
	{
		new(&m_elements[i]) T(m_elements[i+1]);
		m_elements[i+1].~T();
	}
	m_size--;
}

template <typename T>
void array<T>::insert(const T& x, const array_iterator<T>& it)
{
	/*use another iterator to find element indicated by 'it', copy x into this position and shift others*/
	int i;
	int pos;
	array_iterator<T> ti;
	ti=begin();
	while(ti!=end())
	{
		if(ti == it)
			{break;}
		pos++;
		ti++;
	}
	if(m_size<m_reserved_size)
	{
		for(i=m_size;i>pos;i--)
		{
			m_elements[i]=m_elements[i-1];
		}
	new(&m_elements[pos]) T(x);
	m_size++;
	m_reserved_size=m_size;
	}
	else
	{
		T* temp=(T*) malloc(sizeof(T)*m_size+1);
		for(i=0;i<pos;i++)
		{
			new(&temp[i]) T(m_elements[i]);
			m_elements[i].~T();
		}
		new(&temp[pos]) T(x);
		for(i=pos+1;i<m_size;i++)
		{
			new(&temp[i]) T(temp[i-1]);
			m_elements[i-1].~T();
		}
		free(m_elements);
		m_elements=temp;
		m_size++;
		m_reserved_size=m_size;
	}
}

/*array_iterator just points to the current element --- m_current (type T)*/
template <typename T>
array_iterator<T>::array_iterator()
{
	m_current=nullptr;
}

template <typename T>
array_iterator<T>::array_iterator(T* x)	/*passing the member - so directly assign it*/
{
	m_current=x;
}

template <typename T>
array_iterator<T>::array_iterator(const array_iterator& it)	/*copy corresponding member*/
{
	m_current=it.m_current;
}

template <typename T>
T& array_iterator<T>::operator*() const
{
	return *m_current;
}

template <typename T>
array_iterator<T> array_iterator<T>::operator++()
{
	array_iterator<T> temp;
	m_current++;
	temp.m_current=m_current;
	return temp;
}

template <typename T>
array_iterator<T> array_iterator<T>::operator++(int __unused)
{
	array_iterator<T> temp;
	temp.m_current=m_current;
	m_current++;
	return temp;
}

template <typename T>
bool array_iterator<T>::operator != (const array_iterator& it) const
{
	if(m_current!=it.m_current)
		{return true;}
	else
		{return false;}
}

template <typename T>
bool array_iterator<T>::operator == (const array_iterator& it) const
{
	if(m_current==it.m_current)
		{return true;}
	else
		{return false;}
}

//call the destructor on each element and then free
//use in-place in the context of malloc and free
//m_elements=new T[m_reserved_size] and then copy element by element
//malloc - sizeof times n