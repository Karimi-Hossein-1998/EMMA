#pragma once
#include "../Matrix-Decl.hpp"
#include "../LTMatrix-Decl.hpp"

template <Number T>
inline Matrix<T> LTMatrix<T>::ToDense() const
{
	Matrix<T> result(nDim,nDim); if (nDim==0 || data.empty()) return result;
	for (std::uint64_t r{}; r<nDim; ++r) for (std::uint64_t c{}; c<=r; ++c) result[r,c] = data[((r*(r+1))>>1)+c];
	return result;
}
