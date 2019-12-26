#include <TimeMarchingAlgorithmsFactory.hpp>

namespace RedMA
{

std::shared_ptr<TimeMarchingAlgorithm >
TimeMarchingAlgorithmsFactory(const GetPot& datafile,
                              GlobalAssembler* assembler,
                              std::shared_ptr<Epetra_Comm> comm,
                              bool verbose)
{
    bool steady = datafile("time_discretization/steady", false);
    if (steady)
    {
        typedef SteadySolver  ReturnType;
        std::shared_ptr<ReturnType>
                returnPtr(new SteadySolver(datafile,assembler,comm,verbose));

        return returnPtr;
    }

    std::string marchingAlgorithmString =
        datafile("time_discretization/algorithm", "rosenbrock");

    if (!std::strcmp(marchingAlgorithmString.c_str(), "rosenbrock"))
    {
        typedef RosenbrockAlgorithm  ReturnType;
        std::shared_ptr<ReturnType>
                returnPtr(new RosenbrockAlgorithm(datafile,assembler,comm,verbose));

        return returnPtr;
    }
    else if (!std::strcmp(marchingAlgorithmString.c_str(), "backward_euler"))
    {
        typedef BackwardEuler  ReturnType;
        std::shared_ptr<ReturnType>
                returnPtr(new BackwardEuler(datafile,assembler,comm,verbose));

        return returnPtr;
    }
    else
    {
        std::string errorMsg = "Time marching algorithm of type " +
                    marchingAlgorithmString + " is not implemented!";

        throw Exception(errorMsg);
    }
    return nullptr;
}

}  // namespace RedMA
