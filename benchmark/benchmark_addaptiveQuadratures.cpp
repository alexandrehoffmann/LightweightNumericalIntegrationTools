#include <benchmark/benchmark.h>

#include <LNIT/AdaptiveQuadratures.hpp>

static void benchmark_GaussLegendreAdaptiveQuadrature(benchmark::State& state) 
{
	const int k = state.range(0);
	
	const double mu = 0.5;
    const double sigma = 1;

    const auto f = [mu, sigma, k](const double x) -> double
    {
        return (1 / (sigma*std::sqrt(2*M_PI)))*std::exp(-(x-mu)*(x-mu) / (2*sigma*sigma))*std::pow(x, k);
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
	
	const double mu = 0.5;
    const double sigma = 1;

    const auto f = [mu, sigma, k](const double x) -> double
    {
        return (1 / (sigma*std::sqrt(2*M_PI)))*std::exp(-(x-mu)*(x-mu) / (2*sigma*sigma))*std::pow(x, k);
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
	
	const double mu = 0.5;
    const double sigma = 1;

    const auto f = [mu, sigma, k](const double x) -> double
    {
        return (1 / (sigma*std::sqrt(2*M_PI)))*std::exp(-(x-mu)*(x-mu) / (2*sigma*sigma))*std::pow(x, k);
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

BENCHMARK(benchmark_GaussLegendreAdaptiveQuadrature)->DenseRange(0, 80, 2);
BENCHMARK(benchmark_ClensawCurtisAdaptiveQuadrature)->DenseRange(0, 80, 2);
BENCHMARK(benchmark_ClensawCurtisHybridAdaptiveQuadrature)->DenseRange(0, 80, 2);

static void benchmark_RemapGaussLegendreAdaptiveQuadrature(benchmark::State& state) 
{
	const int k = state.range(0);
	
	const double mu = 0.5;
    const double sigma = 1;

    const auto f = [mu, sigma, k](const double x) -> double
    {
        return (1 / (sigma*std::sqrt(2*M_PI)))*std::exp(-(x-mu)*(x-mu) / (2*sigma*sigma))*std::pow(x, k);
    };
    
	LNIT::GaussLegendreAdaptiveQuadrature<double> quad;
	quad.setTol(1.e-11);
	
    for (auto _ : state) 
    {
		double integral = quad.remapAndIntegrate(f);
        benchmark::DoNotOptimize(integral);
    }
    state.counters.insert({"Nit", quad.getNits()});
}

static void benchmark_RemapClensawCurtisAdaptiveQuadrature(benchmark::State& state) 
{
	const int k = state.range(0);
	
	const double mu = 0.5;
    const double sigma = 1;

    const auto f = [mu, sigma, k](const double x) -> double
    {
        return (1 / (sigma*std::sqrt(2*M_PI)))*std::exp(-(x-mu)*(x-mu) / (2*sigma*sigma))*std::pow(x, k);
    };

	LNIT::ClensawCurtisAdaptiveQuadrature<double> quad;
	quad.setTol(1.e-11);
	
    for (auto _ : state) 
    {
		double integral = quad.remapAndIntegrate(f);
        benchmark::DoNotOptimize(integral);
    }
    state.counters.insert({"Nit", quad.getNits()});
}

static void benchmark_RemapClensawCurtisHybridAdaptiveQuadrature(benchmark::State& state) 
{
	const int k = state.range(0);
	
	const double mu = 0.5;
    const double sigma = 1;

    const auto f = [mu, sigma, k](const double x) -> double
    {
        return (1 / (sigma*std::sqrt(2*M_PI)))*std::exp(-(x-mu)*(x-mu) / (2*sigma*sigma))*std::pow(x, k);
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

BENCHMARK(benchmark_RemapGaussLegendreAdaptiveQuadrature)->DenseRange(0, 80, 2);
BENCHMARK(benchmark_RemapClensawCurtisAdaptiveQuadrature)->DenseRange(0, 80, 2);
BENCHMARK(benchmark_RemapClensawCurtisHybridAdaptiveQuadrature)->DenseRange(0, 80, 2);

BENCHMARK_MAIN();
