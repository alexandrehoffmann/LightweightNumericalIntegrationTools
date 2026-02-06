#include <LNIT/AdaptiveQuadratures/GLCCAdaptiveQuadrature.hpp>

namespace LNIT
{

template class GLCCAdaptiveQuadrature<float, float>;
template class GLCCAdaptiveQuadrature<float, double>;
template class GLCCAdaptiveQuadrature<float, long double>;

template class GLCCAdaptiveQuadrature<double, double>;
template class GLCCAdaptiveQuadrature<double, long double>;

template class GLCCAdaptiveQuadrature<long double, long double>;

} // namespace LNIT
