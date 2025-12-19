#include <LNIT/AdaptiveQuadratures/ClenshawCurtisAdaptiveQuadrature.hpp>

namespace LNIT
{

template class ClenshawCurtisAdaptiveQuadrature<double, double>;
template class ClenshawCurtisAdaptiveQuadrature<double, long double>;

template class ClenshawCurtisAdaptiveQuadrature<long double, long double>;

} // namespace LNIT
