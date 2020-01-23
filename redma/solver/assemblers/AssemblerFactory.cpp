#include "AssemblerFactory.hpp"

namespace RedMA
{

SHP(aAssembler)
AssemblerFactory(const GetPot& datafile)
{
    SHP(aAssembler) ret;
    std::string assemblerString = datafile("assemblertype","stokes");

    if (!std::strcmp(assemblerString.c_str(),"stokes"))
        ret.reset(new StokesAssembler(datafile));
    else
        throw new Exception("Time Marching Algorithm is not implemented!");

    return ret;
}

}
