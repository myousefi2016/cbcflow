
#include <string>
#include <iostream>

// TODO: Make a C++ checklib with macros? Other approach?
void check_fail()
{
    std::cerr << "Check failed" << std::endl;
}

void check(bool truth)
{
    if (!truth)
    {
        check_fail();
    }
}

template<typename L, typename R>
void check_eq(L lhs, R rhs);

template<typename T>
void check_eq(const T & lhs, const T & rhs)
{
    if (!(lhs == rhs))
    {
        check_fail();
    }
}

template<>
void check_eq(const std::string & lhs, const char * rhs)
{
    if (!(lhs == std::string(rhs)))
    {
        check_fail();
    }
}

template<typename T>
void check_lt(const T & lhs, const T & rhs)
{
    if (!(lhs < rhs))
    {
        check_fail();
    }
}

template<typename T>
void check_gt(const T & lhs, const T & rhs)
{
    if (!(lhs > rhs))
    {
        check_fail();
    }
}
