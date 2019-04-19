#include "PrintLog.hpp"

namespace RedMA
{

void
printlog(Color outColor, std::string text, bool verbose)
{
    std::string sp = "  ";
    unsigned int hierarchy = 0;
    if (outColor == MAGENTA)
        hierarchy = 0;
    if (outColor == GREEN)
        hierarchy = 1;
    if (outColor == YELLOW)
        hierarchy = 2;

    if(verbose)
    {
        std::cout << "\033[;" << outColor+30 << "m";
        for (int i = 0; i < hierarchy; i++)
            std::cout << sp;
        std::cout << text;
        std::cout << "\033[0m";
    }
}

template < typename T >
std::string
to_string( const T& n )
{
    std::ostringstream stm ;
    stm << n ;
    return stm.str() ;
}

void
printlog(Color outColor, int num, bool verbose)
{
    printlog(outColor, to_string(num), verbose);
}

void
CoutRedirecter::
redirect()
{
    M_prevBuf = std::cout.rdbuf();
    std::cout.rdbuf(M_strCout.rdbuf());
}

std::string
CoutRedirecter::
restore()
{
    std::cout.rdbuf(M_prevBuf);
    return M_strCout.str();
}

}  // PrintLog
