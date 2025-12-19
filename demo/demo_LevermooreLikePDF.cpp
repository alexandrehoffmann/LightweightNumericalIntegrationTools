#include <LNIT/AdaptiveQuadratures.hpp>

int main()
{
	using Scalar = double;
	using Size   = unsigned int;
	
	const auto f = [](const Scalar x) -> Scalar
	{
		//~ return std::exp(1 -x*x + 2*x*x*x - 0.2*x*x*x*x + 0.5*x*x*x*x*x - 0.1*x*x*x*x*x*x);
		return std::exp(5 - 10*x - 5*x*x - 0.01*x*x*x - 0.02*x*x*x*x);
	};
	
	LNIT::ClenshawCurtisHybridAdaptiveQuadrature<Scalar> quad;
	quad.setTol(1.e-11);
	
	fmt::print("Hessian of Levermoor like PDF:\n");
	for (Size i=0; i!=5; ++i)
	{
		for (Size j=0; j!=5; ++j)
		{
			const auto Hfij = [&f, i, j](const Scalar x) -> Scalar
			{
				return f(x)*std::pow(x, i)*std::pow(x,j);
			};
			
			fmt::print("{} ", quad.integrate(Hfij));
		}
		fmt::print("\n");
	}
	
	return EXIT_SUCCESS;
}
