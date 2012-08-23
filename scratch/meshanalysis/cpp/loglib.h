
#include <iostream>

template<typename T>
void log(T text)
{
    std::cout << text << std::endl;
}

template<typename T0, typename T1>
void log2(T0 text0, T1 text1)
{
    std::cout << text0 << "; " << text1 << std::endl;
}

template<typename T0, typename T1, typename T2>
void log3(T0 text0, T1 text1, T2 text2)
{
    std::cout << text0 << "; " << text1 << "; " << text2 << std::endl;
}

template<typename T0, typename T1, typename T2, typename T3>
void log4(T0 text0, T1 text1, T2 text2, T3 text3)
{
    std::cout << text0 << "; " << text1 << "; " << text2 <<  "; " << text3 << std::endl;
}

template<typename T0, typename T1, typename T2, typename T3, typename T4>
void log5(T0 text0, T1 text1, T2 text2, T3 text3, T4 text4)
{
    std::cout << text0 << "; " << text1 << "; " << text2 <<  "; " << text3 <<  "; " << text4 << std::endl;
}
