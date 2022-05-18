//
// Created by betti on 4/28/22.
//

#include "weightFunction.hpp"

namespace RedMA
{
    weightFunction::
    weightFunction()
    {
        setCenter();
        setRadius();
    }

    return_Type
    weightFunction::
    operator()(const Vector3D& position)
    {
        return ((position-M_center).norm()/M_radius < 1) * (0.5 + 0.5 * std::sin(std::atan(1) * 4 * ((position-M_center).norm()/M_radius + 0.5)));
    }

    void
    weightFunction::
    setCenter()
    {
        M_center[0] = -7.99255;
        M_center[1] = 1.97546;
        M_center[2] = 41.3594;
    }

    void
    weightFunction::
    setRadius()
    {
        M_radius = 0.05;
    }
}