#pragma once
#include "../Matrix-Decl.hpp"
#include "../UTMatrix-Decl.hpp"

template <Number T>
inline Matrix<T> UTMatrix<T>::ToDense() const
{
	Matrix<T> result(nDim,nDim); if (nDim==0 || data.empty()) return result;
	for (std::uint64_t r{}; r<nDim; ++r) for (std::uint64_t c{r}; c<nDim; ++c) result[r,c] = data[((r*((nDim<<1)-r+1))>>1)+nDim-c-1];
	return result;
}
