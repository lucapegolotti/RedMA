#include "Chrono.hpp"

namespace RedMA
{

void
Chrono::
start()
{
    M_initialTime = std::chrono::high_resolution_clock::now();
}

double
Chrono::
diff()
{
    using namespace std::chrono;

    high_resolution_clock::time_point curTime = high_resolution_clock::now();

    duration<double> timeSpan = duration_cast<duration<double>>(curTime - M_initialTime);

    return timeSpan.count();
}

}  // namespace RedMA
