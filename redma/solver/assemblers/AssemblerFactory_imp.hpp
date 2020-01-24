#include "AssemblerFactory.hpp"

namespace RedMA
{

template<class InVectorType, class InMatrixType>
SHP(aAssembler<InVectorType AND InMatrixType>)
AssemblerFactory(const GetPot& datafile)
{
    std::shared_ptr<aAssembler<InVectorType, InMatrixType>> ret;
    std::string assemblerString = datafile("assembler/type","stokes");

    if (!std::strcmp(assemblerString.c_str(),"stokes"))
        ret.reset(new StokesAssembler<InVectorType,InMatrixType>(datafile));
    else
        throw new Exception("Time Marching Algorithm is not implemented!");

    return ret;
}

}
