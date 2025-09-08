#include <benchmark/benchmark.h>

#include <LNIT/AdaptiveQuadratures.hpp>

static void benchmark_GaussLegendreAdaptiveQuadrature(benchmark::State& state) 
{
	const int k = state.range(0);
	
	const auto f = [k](const double& x) -> double
	{
		return std::exp(-x*x)*std::pow(x, k);
	};
	
	LNIT::GaussLegendreAdaptiveQuadrature<double> quad;
	quad.setTol(1.e-11);
	
    for (auto _ : state) 
    {
		double integral = quad.integrate(f);
        benchmark::DoNotOptimize(integral);
    }
    state.counters.insert({"Nit", quad.getNits()});
}

static void benchmark_ClensawCurtisAdaptiveQuadrature(benchmark::State& state) 
{
	const int k = state.range(0);
	
	const auto f = [k](const double& x) -> double
	{
		return std::exp(-x*x)*std::pow(x, k);
	};
	
	LNIT::ClensawCurtisAdaptiveQuadrature<double> quad;
	quad.setTol(1.e-11);
	
    for (auto _ : state) 
    {
		double integral = quad.integrate(f);
        benchmark::DoNotOptimize(integral);
    }
    state.counters.insert({"Nit", quad.getNits()});
}

static void benchmark_ClensawCurtisHybridAdaptiveQuadrature(benchmark::State& state) 
{
	const int k = state.range(0);
	
	const auto f = [k](const double& x) -> double
	{
		return std::exp(-x*x)*std::pow(x, k);
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


BENCHMARK(benchmark_GaussLegendreAdaptiveQuadrature)->DenseRange(0, 100, 2);
BENCHMARK(benchmark_ClensawCurtisAdaptiveQuadrature)->DenseRange(0, 100, 2);
BENCHMARK(benchmark_ClensawCurtisHybridAdaptiveQuadrature)->DenseRange(0, 100, 2);

BENCHMARK_MAIN();
