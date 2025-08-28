#include <LNIT/AdaptiveQuadratures.hpp>

inline constexpr double normalDistribRawMoments(const double mu, const double sigma, const unsigned int k)
{
    if (k == 0) { return 1.; }
    if (k == 1) { return mu; }
    return mu*normalDistribRawMoments(mu, sigma, k-1) + (k-1)*sigma*sigma*normalDistribRawMoments(mu, sigma, k-2);
}

int main()
{
	using Scalar = long double;
	using Size   = unsigned int;
	
    const Scalar mu = 0.5;
    const Scalar sigma = 1;

    const auto f = [mu, sigma](const Scalar x) -> Scalar
    {
        return (1 / (sigma*std::sqrt(2*M_PI)))*std::exp(-(x-mu)*(x-mu) / (2*sigma*sigma));
    };

    LNIT::ClensawCurtisAdaptiveQuadrature<Scalar> quad;
    quad.setTol(1.e-10);
    
	std::FILE* fout = stdout;
	fmt::print(fout, "#Estimations of the moments of the normal distribution with mu = {} and sigma = {}\n", mu, sigma);
	fmt::print(fout, "#order has_converged n_iterations estimated analytical absolute_error relative_error\n");

    for (Size k=0; k!=41; ++k)
    {		
		std::FILE* quadOut = std::fopen(fmt::format("Clensaw_Curtis_for_moment_{}.log", k).c_str(), "w");
		quad.setOutput(quadOut);
		
        const auto mf_k = [&f, k](const Scalar x) -> Scalar
        {
            return f(x)*std::pow(x, k);
        };
        const Scalar estimated  = quad.integrate(mf_k);
        const Scalar analytical = normalDistribRawMoments(mu, sigma, k);
        fmt::print(fout, "{} {} {} {} {} {} {}\n", k, quad.hasConverged(), quad.getNits(), estimated, analytical, std::abs(estimated - analytical), std::abs(estimated - analytical) / analytical);
        std::fflush(fout);
        
        std::fclose(quadOut);
    }
    std::fclose(fout);
	
	return EXIT_SUCCESS;
}
