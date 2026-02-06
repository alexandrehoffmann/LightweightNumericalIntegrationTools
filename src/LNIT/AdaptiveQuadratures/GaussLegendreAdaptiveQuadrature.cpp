#include <LNIT/AdaptiveQuadratures/GaussLegendreAdaptiveQuadrature.hpp>

namespace LNIT
{

template class GaussLegendreAdaptiveQuadrature<float, float>;
template class GaussLegendreAdaptiveQuadrature<float, double>;
template class GaussLegendreAdaptiveQuadrature<float, long double>;

template class GaussLegendreAdaptiveQuadrature<double, double>;
template class GaussLegendreAdaptiveQuadrature<double, long double>;

template class GaussLegendreAdaptiveQuadrature<long double, long double>;

} // namespace LNIT
