#include <LNIT/GaussHermiteQuadrature.hpp>

namespace LNIT
{

template class GaussHermiteQuadrature<float, float>;
template class GaussHermiteQuadrature<float, double>;
template class GaussHermiteQuadrature<float, long double>;

template class GaussHermiteQuadrature<double, double>;
template class GaussHermiteQuadrature<double, long double>;

template class GaussHermiteQuadrature<long double, long double>;

} // namespace LNIT
