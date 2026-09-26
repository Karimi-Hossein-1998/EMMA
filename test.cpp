#include "src/MM/typedefs/header.hpp"
#include "src/MM/utility/write.hpp"
#include "src/MM/solvers/ODE/multistep/ab-solver.hpp"
#include <cstddef>
#include <random>

MathEngine::MyFunc GenerateRHS(double alpha)
{
	return [alpha](double t, const MathEngine::dVec& y, MathEngine::dVec& dydt)
	{
		for (size_t i{0}; i<y.size(); ++i)
		{
			dydt[i] = - alpha*y[i];
		}
	};
}

int main(int argc, char* argv[])
{
	double alpha = 0.5;
	if (argc>1) alpha = std::stod(argv[1]);

	size_t n = 10;
	if (argc>2) n = std::stoul(argv[2]);

	std::vector<double> y0; y0.reserve(n);
	std::mt19937_64 rng(41);
	std::uniform_real_distribution<double> urd;

	for (size_t i{0}; i<n; ++i) y0.push_back(urd(rng));

	MathEngine::ODESolverParameters params{
		.derivative        = GenerateRHS(alpha),
		.onStep            = nullptr,
		.initialConditions = y0,
		.t0                = 0.0,
		.t1                = 100.0,
		.dt                = 0.01
	};

	MathEngine::Matrix<double> sol(MathEngine::adams_bashforth_solver(params).solution);

	MathEngine::IO::WriteOptions wo{
		.path      = "test-dir/inner-dir/innermost-dir/the-last-dir/test-sol-"+std::to_string(n)+".csv",
		.separator = ",",
		.header    = {},
		.comment   = {},
		.footer    = {},
		.colWidth  = 20,
		.precision = 15,
		.format    = MathEngine::IO::FPFormat::Scientific,
		.alignment = MathEngine::IO::Alignment::Center,
		.append    = false,
		.binary    = false
	};
	MathEngine::IO::WriteOptions wobin{
		.path      = "test-dir/inner-dir/innermost-dir/the-last-dir/test-sol-"+std::to_string(n)+".bin",
		.separator = ",",
		.header    = {},
		.comment   = {},
		.footer    = {},
		.colWidth  = 20,
		.precision = 15,
		.format    = MathEngine::IO::FPFormat::Scientific,
		.alignment = MathEngine::IO::Alignment::Center,
		.append    = false,
		.binary    = true
	};

	MathEngine::IO::WriteMatrix(sol,wo);
	MathEngine::IO::WriteMatrix(sol,wobin);
	wo.path = "test-dir/inner-dir/innermost-dir/the-last-dir/test-y0-"+std::to_string(n)+".csv"; 
	MathEngine::IO::WriteVector(std::span<const double>(y0), wo);
}




