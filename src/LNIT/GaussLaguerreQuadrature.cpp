#include <LNIT/GaussLaguerreQuadrature.hpp>

namespace LNIT
{

template class GaussLaguerreQuadrature<float, float>;
template class GaussLaguerreQuadrature<float, double>;
template class GaussLaguerreQuadrature<float, long double>;

template class GaussLaguerreQuadrature<double, double>;
template class GaussLaguerreQuadrature<double, long double>;

template class GaussLaguerreQuadrature<long double, long double>;

} // namespace LNIT
