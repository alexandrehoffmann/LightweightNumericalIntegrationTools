#include <LNIT/GaussHermiteQuadrature.hpp>
#include <LNIT/GaussLaguerreQuadrature.hpp>

#include <cmath>
#include <map>

using ArgsTuple = std::tuple<double, double, unsigned int>;

std::map<ArgsTuple, double> memory_normalDistribRawMoments{};

inline double normalDistribRawMoments(const double mu, const double sigma, const unsigned int k)
{	
	using Iterator = typename std::map<ArgsTuple, double>::const_iterator;
	const ArgsTuple argsTuple = std::make_tuple(mu, sigma, k);
	const Iterator  resIt  = memory_normalDistribRawMoments.find(argsTuple);
	
	if (resIt == memory_normalDistribRawMoments.cend())
	{
		if (k == 0) { return 1.; }
		if (k == 1) { return mu; }
		const double res = mu*normalDistribRawMoments(mu, sigma, k-1) + (k-1)*sigma*sigma*normalDistribRawMoments(mu, sigma, k-2);
		memory_normalDistribRawMoments[argsTuple] = res;
		return res;
	}
	else
	{
		return resIt->second;
	}
}

int main()
{
	using Scalar = double;
	using Size   = unsigned int;
	
    const Scalar mu = 0.5;
    const Scalar sigma = 1;

    const auto f = [mu, sigma](const Scalar x) -> Scalar
    {
        return (1 / (sigma*std::sqrt(2*M_PI)))*std::exp(-(x-mu)*(x-mu) / (2*sigma*sigma));
    };

    LNIT::GaussHermiteQuadrature<Scalar>  quadHermite;
    LNIT::GaussLaguerreQuadrature<Scalar> quadLaguerre;
    
	std::FILE* fout = std::fopen("hermite_and_laguerre.log", "w");
	fmt::print(fout, "#Estimations of the moments of the normal distribution with mu = {} and sigma = {}\n", mu, sigma);
	fmt::print(fout, "#order Hermite_estimated Laguerre_estimate analytical Hermite_absolute_error Laguerre_absolute_error Hermite_relative_error Laguerre_relative_error\n");

    for (Size k=0; k!=101; ++k)
    {		
		const auto mf_k = [&f, k](const Scalar x) -> Scalar
        {
            return f(x)*std::pow(x, k);
        };
		
        const Scalar hermiteEstimated  = quadHermite.integrate(mf_k);
        const Scalar laguerreEstimated = quadLaguerre.integrateLeftInfinite(mf_k) + quadLaguerre.integrateRightInfinite(mf_k);
        const Scalar analytical = normalDistribRawMoments(mu, sigma, k);
        fmt::print(fout, "{} {} {} {} {} {} {} {}\n", k, hermiteEstimated, laguerreEstimated, analytical, std::abs(hermiteEstimated - analytical), std::abs(laguerreEstimated - analytical), std::abs(hermiteEstimated - analytical) / analytical, std::abs(laguerreEstimated - analytical) / analytical);
        std::fflush(fout);
    }
    std::fclose(fout);
	
	return EXIT_SUCCESS;
}
