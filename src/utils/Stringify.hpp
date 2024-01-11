#pragma once

#include <iomanip>
#include <iostream>
#include <string>

namespace utils
{

template <class T>
class StringifyImpl
{
public:
    template <class U>
    friend std::ostream& operator<<(std::ostream& os, const StringifyImpl<U>& prop);

public:
    StringifyImpl(const std::string& name, const T& value);

private:
    std::string m_name;
    T m_value;
};

template <class T>
inline
std::ostream& operator<<(std::ostream& os, const StringifyImpl<T>& prop)
{
    os << std::left << std::setw(40) << std::setfill('.') << prop.m_name;
    os << std::right << std::setw(40) << std::setfill('.') << prop.m_value;
    return os;
}

template <class T>
inline
StringifyImpl<T>::StringifyImpl(const std::string& name, const T& value)
    : m_name (name)
    , m_value (value)
{
}

template <class T>
inline
static StringifyImpl<T> stringify(const std::string& name, const T& value)
{
    return StringifyImpl<T> (name, value);
}

}  // utils
