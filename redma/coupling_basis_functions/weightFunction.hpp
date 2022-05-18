//
// Created by betti on 4/28/22.
//

#include <redma/RedMA.hpp>

#ifndef WEIGHTFUNCTION_HPP
#define WEIGHTFUNCTION_HPP

namespace RedMA
{
    typedef double return_Type;
    typedef LifeV::VectorSmall<3> Vector3D;

class weightFunction
    {
    public:
        typedef double return_Type;
        weightFunction();
        return_Type operator()(const Vector3D& position);

        void setRadius();
        void setCenter();

        return_Type getRadius() { return M_radius; };
        Vector3D getCenter() { return M_center; };
    private:
    return_Type M_radius;
    Vector3D M_center;
    };
}

#endif //WEIGHTFUNCTION_HPP
