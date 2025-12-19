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

static void benchmark_ClenshawCurtisAdaptiveQuadrature(benchmark::State& state) 
{
	const int k = state.range(0);
	
	const double mu = 0.5;
    const double sigma = 1;

    const auto f = [mu, sigma, k](const double x) -> double
    {
        return (1 / (sigma*std::sqrt(2*M_PI)))*std::exp(-(x-mu)*(x-mu) / (2*sigma*sigma))*std::pow(x, k);
    };
	
	LNIT::ClenshawCurtisAdaptiveQuadrature<double> quad;
	quad.setTol(1.e-11);
	
    for (auto _ : state) 
    {
		double integral = quad.integrate(f);
        benchmark::DoNotOptimize(integral);
    }
    state.counters.insert({"Nit", quad.getNits()});
}

static void benchmark_ClenshawCurtisHybridAdaptiveQuadrature(benchmark::State& state) 
{
	const int k = state.range(0);
	
	const double mu = 0.5;
    const double sigma = 1;

    const auto f = [mu, sigma, k](const double x) -> double
    {
        return (1 / (sigma*std::sqrt(2*M_PI)))*std::exp(-(x-mu)*(x-mu) / (2*sigma*sigma))*std::pow(x, k);
    };
	
	LNIT::ClenshawCurtisHybridAdaptiveQuadrature<double> quad;
	quad.setTol(1.e-11);
	
    for (auto _ : state) 
    {
		double integral = quad.integrate(f);
        benchmark::DoNotOptimize(integral);
    }
    state.counters.insert({"Nit", quad.getNits()});
}

static void benchmark_GLCCAdaptiveQuadrature(benchmark::State& state) 
{
	const int k = state.range(0);
	
	const double mu = 0.5;
    const double sigma = 1;

    const auto f = [mu, sigma, k](const double x) -> double
    {
        return (1 / (sigma*std::sqrt(2*M_PI)))*std::exp(-(x-mu)*(x-mu) / (2*sigma*sigma))*std::pow(x, k);
    };
	
	LNIT::GLCCAdaptiveQuadrature<double> quad;
	quad.setTol(1.e-11);
	
    for (auto _ : state) 
    {
		double integral = quad.integrate(f);
        benchmark::DoNotOptimize(integral);
    }
    state.counters.insert({"Nit", quad.getNits()});
}

BENCHMARK(benchmark_GaussLegendreAdaptiveQuadrature)->DenseRange(0, 100, 2);
BENCHMARK(benchmark_ClenshawCurtisAdaptiveQuadrature)->DenseRange(0, 100, 2);
BENCHMARK(benchmark_ClenshawCurtisHybridAdaptiveQuadrature)->DenseRange(0, 100, 2);
BENCHMARK(benchmark_GLCCAdaptiveQuadrature)->DenseRange(0, 100, 2);

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

static void benchmark_RemapClenshawCurtisAdaptiveQuadrature(benchmark::State& state) 
{
	const int k = state.range(0);
	
	const double mu = 0.5;
    const double sigma = 1;

    const auto f = [mu, sigma, k](const double x) -> double
    {
        return (1 / (sigma*std::sqrt(2*M_PI)))*std::exp(-(x-mu)*(x-mu) / (2*sigma*sigma))*std::pow(x, k);
    };

	LNIT::ClenshawCurtisAdaptiveQuadrature<double> quad;
	quad.setTol(1.e-11);
	
    for (auto _ : state) 
    {
		double integral = quad.remapAndIntegrate(f);
        benchmark::DoNotOptimize(integral);
    }
    state.counters.insert({"Nit", quad.getNits()});
}

static void benchmark_RemapClenshawCurtisHybridAdaptiveQuadrature(benchmark::State& state) 
{
	const int k = state.range(0);
	
	const double mu = 0.5;
    const double sigma = 1;

    const auto f = [mu, sigma, k](const double x) -> double
    {
        return (1 / (sigma*std::sqrt(2*M_PI)))*std::exp(-(x-mu)*(x-mu) / (2*sigma*sigma))*std::pow(x, k);
    };
	
	LNIT::ClenshawCurtisHybridAdaptiveQuadrature<double> quad;
	quad.setTol(1.e-11);
	
    for (auto _ : state) 
    {
		double integral = quad.remapAndIntegrate(f);
        benchmark::DoNotOptimize(integral);
    }
    state.counters.insert({"Nit", quad.getNits()});
}

static void benchmark_RemapGLCCAdaptiveQuadrature(benchmark::State& state) 
{
	const int k = state.range(0);
	
	const double mu = 0.5;
    const double sigma = 1;

    const auto f = [mu, sigma, k](const double x) -> double
    {
        return (1 / (sigma*std::sqrt(2*M_PI)))*std::exp(-(x-mu)*(x-mu) / (2*sigma*sigma))*std::pow(x, k);
    };
	
	LNIT::GLCCAdaptiveQuadrature<double> quad;
	quad.setTol(1.e-11);
	
    for (auto _ : state) 
    {
		double integral = quad.remapAndIntegrate(f);
        benchmark::DoNotOptimize(integral);
    }
    state.counters.insert({"Nit", quad.getNits()});
}

BENCHMARK(benchmark_RemapGaussLegendreAdaptiveQuadrature)->DenseRange(0, 100, 2);
BENCHMARK(benchmark_RemapClenshawCurtisAdaptiveQuadrature)->DenseRange(0, 100, 2);
BENCHMARK(benchmark_RemapClenshawCurtisHybridAdaptiveQuadrature)->DenseRange(0, 100, 2);
BENCHMARK(benchmark_RemapGLCCAdaptiveQuadrature)->DenseRange(0, 100, 2);

BENCHMARK_MAIN();
