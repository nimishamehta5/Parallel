#pragma once

#include <iostream>

//students should not change this

using size_t = std::size_t;

template<typename T> class array_iterator;

template<typename T>
class array {
public: 
    //default constructor
    array();

    //initialize array with elements in initializer
    array(std::initializer_list<T>);

    //copy constructor
    array(const array&);

    //move constructor
    array(array&&);

    //construct with initial "reserved" size
    array(size_t);

    //construct with n copies of t
    array(size_t n, const T& t);

    //destructor
    ~array();

    //ensure enough memory for n elements - if n is smalller than m_res_size do nothing, else grow the dynamic array size
    void reserve(size_t n);

    //add to end of vector
    void push_back(const T&);

    //add to front of vector
    void push_front(const T&);

    //remove last element
    void pop_back();

    //remove first element
    void pop_front();

    //return reference to first element
    T& front() const;

    //return reference to last element
    T& back() const;

    //return reference to specified element
    const T& operator[](size_t) const;

    //return reference to specified element
    T& operator[](size_t);

    //return number of elements
    size_t length() const;

    //returns true if empty
    bool empty() const;

    //remove all elements
    void clear();

    //obtain iterator to first element
    array_iterator<T> begin() const;

    //obtain iterator to one beyond element
    array_iterator<T> end() const;

    //remove specified element
    void erase(const array_iterator<T>&);

    //insert element right before itr
    void insert(const T&, const array_iterator<T>&);

private:
    T* m_elements;              //points to actual elements
    size_t m_size;              //number of elements
    size_t m_reserved_size;     //number of reserved elements
};

template<typename T>
class array_iterator {
public:
    array_iterator();
    array_iterator(T*);
    array_iterator(const array_iterator&);
    T& operator*() const;
    array_iterator operator++();
    array_iterator operator++(int __unused);
    bool operator != (const array_iterator&) const;
    bool operator == (const array_iterator&) const;

private:
    T* m_current;
    // I want my array class to be able to access m_current even though it's private
    // 'friend' allows this to happen
    friend class array<T>;
};