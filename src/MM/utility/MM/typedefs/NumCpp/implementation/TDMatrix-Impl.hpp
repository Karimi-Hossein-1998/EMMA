#pragma once
#include "../Matrix-Decl.hpp"
#include "../TDMatrix-Decl.hpp"

template <Number T>
inline Matrix<T> TDMatrix<T>::ToDense() const
{
	Matrix<T> result(nDim,nDim); if (nDim==0 || diag.empty()) return result;
	result[0,0] = diag[0]; result[0,1] = up[0];
	for (std::uint64_t r{1}; r<nDim-1; ++r)
	{
		result[r,r] = diag[r]; result[r,r-1] = lo[r-1]; result[r,r+1] = up[r-1];
	}
	result[nDim-1,nDim-2] = lo[nDim-2]; result[nDim-1,nDim-1] = diag[nDim-1];
	return result;
}

template <Number T>
inline void TDSolveInPlace(const TDMatrix<T>& SystemMatrix, Matrix<T>& RHSMatrix)
{
	std::uint64_t NumEqs = SystemMatrix.Rows(); std::uint64_t NumRhs = RHSMatrix.Cols();
	if (NumEqs!=RHSMatrix.Rows() || NumRhs==0) {std::println("Dimension mismatch... Quitting!"); return;}
	if (NumEqs==0) return;
	if (NumEqs==1) {for (std::uint64_t i{}; i<NumRhs; ++i) RHSMatrix[0,i] /= SystemMatrix[0,0]; return;}
	const auto& lo = SystemMatrix.GetLo();
	const auto& diag = SystemMatrix.GetDiag();
	const auto& up = SystemMatrix.GetUp();
	T denom = diag[0]; std::vector<T> c_buffer(NumEqs-1,T{});
	if (denom==T{}) {std::println("Encountered ZERO pivot... Quitting!"); return;}
	c_buffer[0] = up[0]/denom; for (std::uint64_t i{}; i<NumRhs; ++i) RHSMatrix[0,i] /= denom;
	for (std::uint64_t i{1}; i<NumEqs; ++i)
	{
		denom = diag[i] - lo[i-1]*c_buffer[i-1];
		if (denom==T{}) {std::println("Encountered ZERO pivot... Quitting!"); return;}
		if (i<NumEqs-1) c_buffer[i] = up[i] / denom;
		for (std::uint64_t j{}; j<NumRhs; ++j) RHSMatrix[i,j] = (RHSMatrix[i,j]-lo[i-1]*RHSMatrix[i-1,j]) / denom;
	}
	for (std::int64_t i{static_cast<std::int64_t>(NumEqs-2)}; i>=0; --i) for (std::uint64_t j{}; j<NumRhs; ++j) RHSMatrix[i,j] -= c_buffer[i]*RHSMatrix[i+1,j];
}

template <Number T>
inline Matrix<T> TDSolve(const TDMatrix<T>& SystemMatrix, const Matrix<T>& RHSMatrix)
{
	auto b=RHSMatrix; TDSolveInPlace(SystemMatrix,b); return b;
}
