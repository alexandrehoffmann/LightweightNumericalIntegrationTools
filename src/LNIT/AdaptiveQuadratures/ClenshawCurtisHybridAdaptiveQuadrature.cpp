#include <LNIT/AdaptiveQuadratures/ClenshawCurtisHybridAdaptiveQuadrature.hpp>

namespace LNIT
{

template class ClenshawCurtisHybridAdaptiveQuadrature<float, float>;
template class ClenshawCurtisHybridAdaptiveQuadrature<float, double>;
template class ClenshawCurtisHybridAdaptiveQuadrature<float, long double>;

template class ClenshawCurtisHybridAdaptiveQuadrature<double, double>;
template class ClenshawCurtisHybridAdaptiveQuadrature<double, long double>;

template class ClenshawCurtisHybridAdaptiveQuadrature<long double, long double>;

} // namespace LNIT
