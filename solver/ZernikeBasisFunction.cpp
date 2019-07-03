#include <ZernikeBasisFunction.hpp>

namespace RedMA
{

ZernikeBasisFunction::
ZernikeBasisFunction(const GeometricFace& face,
                     unsigned int nMax) :
  BasisFunctionFunctor(face)
{
    M_nMax = nMax;

    double radius = face.M_radius;

    fillFactorials(nMax);

    for (int n = 0; n <= nMax; n++)
    {
        for (int m = -n; m <= n; m++)
        {
            if (!((n - std::abs(m)) % 2))
            {
                M_ms.push_back(m);
                M_ns.push_back(n);

                std::vector<int> curListCoefs;
                unsigned int mult = -1;
                unsigned int num;
                unsigned int den;
                for (int k = 0; k <= (n-std::abs(m))/2; k++)
                {
                    mult = mult * (-1);
                    num = mult * M_factorials[n - k];
                    den = M_factorials[k] * M_factorials[(n+std::abs(m))/2 - k] *
                          M_factorials[(n-std::abs(m))/2-k];
                    curListCoefs.push_back(num / den);
                }
                M_polyCoefs.push_back(curListCoefs);
            }
        }
    }

    M_nBasisFunctions = M_polyCoefs.size();

    Vector3D& normal = M_face.M_normal;

    // arbitrary vector to measure the angle
    // we check just the first component of the normal because we trust that
    // it is unitary (hence if normal[1] == 1 => normal = (1,0,0))
    if (std::abs(std::abs(normal[0]) - 1.0) > 1e-12)
    {
        M_e[0] = 1.0; M_e[1] = 0.0; M_e[2] = 0.0;
    }
    else
    {
        M_e[0] = 0.0; M_e[1] = 1.0; M_e[2] = 0.0;
    }
    // project the vector onto the face and orthonormalize
    M_e = M_e - M_e.dot(normal) * normal;
    M_e = M_e / M_e.norm();
}

void
ZernikeBasisFunction::
setIndex(const unsigned int& index)
{
    BasisFunctionFunctor::setIndex(index);
    if (M_ms[index] < 0)
        M_curFunction = (double(*)(double))&std::sin;
    else
        M_curFunction = (double(*)(double))&std::cos;
}

void
ZernikeBasisFunction::
fillFactorials(unsigned int nMax)
{
    unsigned int curN = 1;
    M_factorials.push_back(curN);
    for (int i = 1; i <= nMax; i++)
    {
        curN *= i;
        M_factorials.push_back(curN);
    }
}

ZernikeBasisFunction::return_Type
ZernikeBasisFunction::
operator()(const Vector3D& pos)
{
    Vector3D& center = M_face.M_center;
    Vector3D& normal = M_face.M_normal;

    Vector3D diff = pos - center;
    double r = diff.norm();

    double theta = std::acos(diff.dot(M_e) / (diff.norm()));
    Vector3D crossProduct = diff.cross(M_e);

    if (crossProduct.dot(normal) < 0)
        theta = theta + M_PI;

    double returnVal = 0;

    unsigned int n = M_ns[M_index];
    unsigned int m = std::abs(M_ms[M_index]);
    std::vector<int>& coefList = M_polyCoefs[M_index];
    for (int k = 0; k <= (n-m)/2; k++)
    {
        returnVal += coefList[k] * std::pow(r,n-2*k);
    }
    returnVal *= M_curFunction(m * theta);
    return returnVal;
}

}  // namespace RedMA
