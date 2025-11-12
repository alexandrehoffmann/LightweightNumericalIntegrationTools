#include <benchmark/benchmark.h>

#include <LNIT/AdaptiveQuadratures.hpp>

static void benchmark_integrateLevermoreLikePDF(benchmark::State& state) 
{
	const auto f = [](const double x) -> double
	{
		return std::exp(1 -x*x + 2*x*x*x - 0.2*x*x*x*x + 0.5*x*x*x*x*x - 0.1*x*x*x*x*x*x);
	};
	
	LNIT::ClensawCurtisHybridAdaptiveQuadrature<double> quad;
	quad.setTol(1.e-11);
	
    for (auto _ : state) 
    {
		double integral = quad.integrate(f);
        benchmark::DoNotOptimize(integral);
    }
    state.counters.insert({"Nit", quad.getNits()});
}

static void benchmark_remapAndIntegrateLevermoreLikePDF(benchmark::State& state) 
{
	const auto f = [](const double x) -> double
	{
		return std::exp(1 -x*x + 2*x*x*x - 0.2*x*x*x*x + 0.5*x*x*x*x*x - 0.1*x*x*x*x*x*x);
	};
	
	LNIT::ClensawCurtisHybridAdaptiveQuadrature<double> quad;
	quad.setTol(1.e-11);
	
    for (auto _ : state) 
    {
		double integral = quad.remapAndIntegrate(f);
        benchmark::DoNotOptimize(integral);
    }
    state.counters.insert({"Nit", quad.getNits()});
}

BENCHMARK(benchmark_integrateLevermoreLikePDF);
BENCHMARK(benchmark_remapAndIntegrateLevermoreLikePDF);

BENCHMARK_MAIN();
