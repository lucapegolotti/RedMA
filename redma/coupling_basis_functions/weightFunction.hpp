//
// Created by betti on 4/28/22.
//

#ifndef WEIGHTFUNCTION_H
#define WEIGHTFUNCTION_H

namespace RedMA {
    class weightFunction {
        typedef LifeV::VectorSmall<3> Vector3D;
    public:
        weightFunction() {}

        double operator()(const Vector3D& position);
        double getRadius() { return M_radius};
        double getCenter() { return M_center};
    private:
        double M_radius;
        Vector3D M_center;
    };
}

#endif //WEIGHTFUNCTION_H
