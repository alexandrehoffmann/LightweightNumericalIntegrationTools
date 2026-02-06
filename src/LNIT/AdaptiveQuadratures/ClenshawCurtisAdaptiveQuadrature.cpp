#include <LNIT/AdaptiveQuadratures/ClenshawCurtisAdaptiveQuadrature.hpp>

namespace LNIT
{

template class ClenshawCurtisAdaptiveQuadrature<float, float>;
template class ClenshawCurtisAdaptiveQuadrature<float, double>;
template class ClenshawCurtisAdaptiveQuadrature<float, long double>;

template class ClenshawCurtisAdaptiveQuadrature<double, double>;
template class ClenshawCurtisAdaptiveQuadrature<double, long double>;

template class ClenshawCurtisAdaptiveQuadrature<long double, long double>;

} // namespace LNIT
