#pragma once
#include "../Matrix-Decl.hpp"
#include "../LTMatrix-Decl.hpp"
#include "../UTMatrix-Decl.hpp"

/***************************************************
 *            _________              ___________  *
 *           / _______ \             |___   ___|  *
 *          / /       \ \                | |      *
 *         / /         \ \               | |      *
 *        / /           \ \              | |      *
 *       | |             | |             | |      *
 *       | |             |_|             | |      *
 *       | |                             | |      *
 *       | |       _________    _        | |      *
 *       | |       |_____  |   | |       | |      *
 *        \ \           / /    | |       | |      *
 *         \ \         / /      \ \     / /       *
 *          \ \_______/ /        \ \___/ /        *
 *           \_________/          \_____/         *
 *                                                *
 ***************************************************/

// Gauss-Jordan Solver (Linear Equation System) (In Place Solver)
template <FPNumber T>
inline void GaussJordanInPlace(Matrix<T>& SystemMatrix, Matrix<T>& RHSMatrix)
{
	const std::uint64_t NumEqs(SystemMatrix.Rows()); const std::uint64_t NumUnKnowns(RHSMatrix.Rows()); const std::uint64_t NumRhs(RHSMatrix.Cols());
	if (NumEqs!=SystemMatrix.Cols()) {std::print("Not a standard problem! Number of Eqs ({}) is not equal to number of Unknowns ({})!\nTherefore quitting...\n",NumEqs,SystemMatrix.Cols()); return;}
	if (NumEqs!=NumUnKnowns) {std::print("Not a standard problem! Number of Unknowns ({}) is not equal to number of Eqs ({})!\nTherefore quitting...\n",NumUnKnowns,NumEqs); return;}

	// Counting pivots
	std::vector<std::uint32_t> is_pivoted(NumEqs,0);
	// Recording pivot locations
	std::vector<std::uint64_t> row_pivot_indices(NumEqs);
	std::vector<std::uint64_t> col_pivot_indices(NumEqs);
	T pivot_reciprocal; T pivot_value{}; T max_value{}; uint64_t pivot_row{}; std::uint64_t pivot_col{}; T factor{};
	for (std::uint64_t step{}; step<NumEqs; ++step)
	{
		// Finding the largest value to pivot (pivot point)
		max_value = 0.0;
		for (std::uint64_t r{}; r<NumEqs; ++r)
		{
			if (is_pivoted[r]!=1)
			{
				for (std::uint64_t c{}; c<NumEqs; ++c)
				{
					if (is_pivoted[c]==0)
					{
						std::float64_t abs_rc = std::abs(SystemMatrix[r][c]);
						if (abs_rc>max_value)
						{
							max_value = abs_rc; pivot_row = r; pivot_col = c;
						}
					}
				}
			}
		}
		++is_pivoted[pivot_col];
		// Normalize the diagonal
		if (pivot_row!=pivot_col)
		{
			SystemMatrix.SwapRows(pivot_row,pivot_col); RHSMatrix.SwapRows(pivot_row,pivot_col);
		}
		row_pivot_indices[step] = pivot_row; col_pivot_indices[step] = pivot_col;
		pivot_value = SystemMatrix[pivot_col,pivot_col];
		if (std::abs(pivot_value)<=std::numeric_limits<T>::epsilon()) {std::println("Impossible to pivot... Quiting gracefully!"); return;}
		pivot_reciprocal = static_cast<T>(1.0/pivot_value);
		for (std::uint64_t c{}; c<NumEqs; ++c) SystemMatrix[pivot_col,c] *= pivot_reciprocal;
		for (std::uint64_t c{}; c<NumRhs; ++c) RHSMatrix[pivot_col,c] *= pivot_reciprocal;
		SystemMatrix[pivot_col][pivot_col] = static_cast<T>(1.0);

		// Reduce the rows to 0
		for (std::uint64_t r{}; r<NumEqs; ++r)
		{
			if (r!=pivot_col)
			{
				factor = SystemMatrix[r,pivot_col];
				for (std::uint64_t c{}; c<NumEqs; ++c) SystemMatrix[r,c] -= SystemMatrix[pivot_col,c]*factor;
				SystemMatrix[r,pivot_col] = T{};
				for (std::uint64_t c{}; c<NumRhs; ++c) RHSMatrix[r,c] -= RHSMatrix[pivot_col,c]*factor;
			}
		}
	}
}

// Gauss-Jordan Solver (Linear Equation System)
template <typename T>
inline Matrix<T> GaussJordan(const Matrix<T>& SystemMatrix, const Matrix<T>& RHSMatrix) {Matrix<T> A(SystemMatrix); Matrix<T> b(RHSMatrix); GaussJordanInPlace(A,b); return b;}

/*************************
 *  _       ___      ___ *
 * | |      | |      | | *
 * | |      | |      | | *
 * | |      | |      | | *
 * | |      | |      | | *
 * | |      | |      | | *
 * | |       \ \    / /  *
 * | |____    \ \__/ /   *
 * |_____|     \____/    *
 *                       *
 *************************/

// L-U Decomposer
// Two Matrix Decompose (Do not recommend using this)
template <FPNumber T>
inline std::pair<Matrix<T>,Matrix<T>> LUDecomposeNormal(const Matrix<T>& SystemMatrix)
{
	std::uint64_t rows = SystemMatrix.Rows(); std::uint64_t cols = SystemMatrix.Cols();
	Matrix<T> L(rows,cols); Matrix<T> U(rows,cols);
	if (rows!=cols) {std::println("Matrix is not square... Gracefully quitting!"); return {L,U};}
	// U[0,0] = SystemMatrix[0,0];
	for (std::uint64_t i{}; i<rows; ++i) {U[0,i]=SystemMatrix[0,i]; L[i,0]=SystemMatrix[i,0]/U[0,0];}
	for (std::uint64_t i{}; i<rows; ++i) L[i,i] = 1.0;
	for (std::uint64_t r{1}; r<rows; ++r)
	{
		for (std::uint64_t c{r}; c<cols; ++c)
		{
			U[r,c] = SystemMatrix[r,c];
			for (std::uint64_t m{}; m<r; ++m)
			{
				U[r,c] -= L[r,m]*U[m,c];
			}
		}
		for (std::uint64_t c{r+1}; c<cols; ++c)
		{
			L[c,r] = SystemMatrix[c,r];
			for (std::uint64_t m{}; m<r; ++m)
			{
				L[c,r] -= L[c,m]*U[m,r];
			}
			L[c,r] /= U[r,r];
		}
	}
	return {L,U};
}
// Single Matrix Decompose
template <FPNumber T>
inline void LUDecomposeInPlace(Matrix<T>& SystemMatrix)
{
	std::uint64_t rows = SystemMatrix.Rows(); std::uint64_t cols = SystemMatrix.Cols();
	if (rows!=cols) {std::println("Matrix is not square... Gracefully quitting!"); return;}
	for (std::uint64_t i{1}; i<rows; ++i) SystemMatrix[i,0]/=SystemMatrix[0,0];
	for (std::uint64_t r{1}; r<rows; ++r)
	{
		for (std::uint64_t c{r}; c<cols; ++c)
		{
			for (std::uint64_t m{}; m<r; ++m)
			{
				SystemMatrix[r,c] -= SystemMatrix[r,m]*SystemMatrix[m,c];
			}
		}
		for (std::uint64_t c{r+1}; c<cols; ++c)
		{
			for (std::uint64_t m{}; m<r; ++m)
			{
				SystemMatrix[c,r] -= SystemMatrix[c,m]*SystemMatrix[m,r];
			}
			SystemMatrix[c,r] /= SystemMatrix[r,r];
		}
	}
}
// Single Matrix Decompose
template <FPNumber T>
inline Matrix<T> LUDecomposeCompact(const Matrix<T>& SystemMatrix)
{
	std::uint64_t rows = SystemMatrix.Rows(); std::uint64_t cols = SystemMatrix.Cols();
	Matrix<T> LU(SystemMatrix);
	if (rows!=cols) {std::println("Matrix is not square... Gracefully quitting!"); return LU;}
	for (std::uint64_t i{1}; i<rows; ++i) LU[i,0]/=LU[0,0];
	for (std::uint64_t r{1}; r<rows; ++r)
	{
		for (std::uint64_t c{r}; c<cols; ++c)
		{
			for (std::uint64_t m{}; m<r; ++m)
			{
				LU[r,c] -= LU[r,m]*LU[m,c];
			}
		}
		for (std::uint64_t c{r+1}; c<cols; ++c)
		{
			for (std::uint64_t m{}; m<r; ++m)
			{
				LU[c,r] -= LU[c,m]*LU[m,r];
			}
			LU[c,r] /= LU[r,r];
		}
	}
	return LU;
}
// L-U Decomposer
template <FPNumber T>
inline std::pair<LTMatrix<T>,UTMatrix<T>> LUDecompose(const Matrix<T>& SystemMatrix)
{
	std::uint64_t rows = SystemMatrix.Rows(); std::uint64_t cols = SystemMatrix.Cols();
	LTMatrix<T> L(rows); UTMatrix<T> U(cols);
	if (rows!=cols) {std::println("Matrix is not square... Gracefully quitting!"); return {L,U};}
	for (std::uint64_t i{}; i<rows; ++i) L[i,i] = 1.0;
	// U[0,0] = SystemMatrix[0,0];
	for (std::uint64_t i{}; i<rows; ++i) {U[0,i]=SystemMatrix[0,i]; L[i,0]=SystemMatrix[i,0]/U[0,0];}
	for (std::uint64_t r{1}; r<rows; ++r)
	{
		for (std::uint64_t c{r}; c<cols; ++c)
		{
			U[r,c] = SystemMatrix[r,c];
			for (std::uint64_t m{}; m<r; ++m)
			{
				U[r,c] -= L[r,m]*U[m,c];
			}
		}
		for (std::uint64_t c{r+1}; c<cols; ++c)
		{
			L[c,r] = SystemMatrix[c,r];
			for (std::uint64_t m{}; m<r; ++m)
			{
				L[c,r] -= L[c,m]*U[m,r];
			}
			L[c,r] /= U[r,r];
		}
	}
	return {L,U};
}

/************************************************************************
 *      ___                  _         _       _                        *
 *     /___\                | |       | |     | |                       *
 *   //     \\              | |       | |     | |                       *
 *  //       \\             | |       | |     | |                       *
 *  ||       ||             | |       | |     | |                       *
 *   \\      ||             | |____    \|\___/|/                        *
 *    \\                    |______|    \_____/                         *
 *      \\             _         _                      _           _   *
 *       \\           /_\       /||    \\     /\       /_\     ||  /_|  *
 *         \\       //   \\     /||    \\    //\\    //   \\   || //    *
 *  ||      \\    //       \\    ||    \\    \\//   |/_____\|  |///     *
 *  ||       ||  ||         ||   ||     \\    XX//  |_______|  ||/      *
 *  \\       //   \\       //    ||      \\  //\/   \\         ||       *
 *   \\_____//      \\___//      ||__     \\//       \\___//   ||       *
 *     \___/          \_/        |__|      \/          \_/     ||       *
 *                                                                      *
 ************************************************************************/

// Solve by L-U Decomposed Matrices (forward + backward substitution) (in-place)
template <FPNumber T>
inline void LUSolveInPlace(const Matrix<T>& L, const Matrix<T>& U, Matrix<T>& RHSMatrix)
{
	std::uint64_t NumEqs = L.Rows(); std::uint64_t NumRhs = RHSMatrix.Cols(); std::uint64_t NumUnKnowns = RHSMatrix.Rows(); std::int64_t NumEqsMinus1 = NumEqs-1;
	if (NumEqs != NumUnKnowns || NumEqs != L.Cols() || NumEqs != U.Rows() || NumEqs != U.Cols()) {std::println("Dimensions mismatched... Quitting!"); return;}
	for (std::uint64_t r{1}; r<NumEqs; ++r)
		for (std::uint64_t c{}; c<r; ++c)
			for (std::uint64_t i{}; i<NumRhs; ++i)
				RHSMatrix[r,i] -= L[r,c]*RHSMatrix[c,i];
	std::vector<T> UrrInv; UrrInv.reserve(NumEqs); for (std::uint64_t i{}; i<NumEqs; ++i) UrrInv.push_back(1.0/U[i,i]);
	for (std::uint64_t i{}; i<NumRhs; ++i)
	{
		RHSMatrix[NumEqsMinus1,i] *= UrrInv[NumEqsMinus1];
	}
	for (std::int64_t r{NumEqsMinus1-1}; r>-1; --r)
	{
		for (std::int64_t c{r+1}; c<NumEqs; ++c)
		{
			for (std::uint64_t i{}; i<NumRhs; ++i)
			{
				RHSMatrix[r,i] -= U[r,c]*RHSMatrix[c,i];
			}
		}
		for (std::uint64_t i{}; i<NumRhs; ++i)
		{
			RHSMatrix[r,i] *= UrrInv[r];
		}
	}
}

// Solve by L-U Decomposed Matrices
template <FPNumber T>
inline Matrix<T> LUSolve(const Matrix<T>& L, const Matrix<T>& U, const Matrix<T>& RHSMatrix) {Matrix<T> b(RHSMatrix); LUSolveInPlace(L,U,b); return b;}

// Solve by L-U Decomposed Matrices (forward + backward substitution) (in-place)
template <FPNumber T>
inline void LUSolveInPlace(const LTMatrix<T>& L, const UTMatrix<T>& U, Matrix<T>& RHSMatrix)
{
	std::uint64_t NumEqs = L.Rows(); std::uint64_t NumRhs = RHSMatrix.Cols(); std::uint64_t NumUnKnowns = RHSMatrix.Rows(); std::int64_t NumEqsMinus1 = NumEqs-1;
	if (NumEqs != NumUnKnowns || NumEqs != L.Cols() || NumEqs != U.Rows() || NumEqs != U.Cols()) {std::println("Dimensions mismatched... Quitting!"); return;}
	for (std::uint64_t r{1}; r<NumEqs; ++r)
		for (std::uint64_t c{}; c<r; ++c)
			for (std::uint64_t i{}; i<NumRhs; ++i)
				RHSMatrix[r,i] -= L[r,c]*RHSMatrix[c,i];
	std::vector<T> UrrInv; UrrInv.reserve(NumEqs); for (std::uint64_t i{}; i<NumEqs; ++i) UrrInv.push_back(1.0/U[i,i]);
	for (std::uint64_t i{}; i<NumRhs; ++i)
	{
		RHSMatrix[NumEqsMinus1,i] *= UrrInv[NumEqsMinus1];
	}
	for (std::int64_t r{NumEqsMinus1-1}; r>-1; --r)
	{
		for (std::int64_t c{r+1}; c<NumEqs; ++c)
		{
			for (std::uint64_t i{}; i<NumRhs; ++i)
			{
				RHSMatrix[r,i] -= U[r,c]*RHSMatrix[c,i];
			}
		}
		for (std::uint64_t i{}; i<NumRhs; ++i)
		{
			RHSMatrix[r,i] *= UrrInv[r];
		}
	}
}

// Solve by L-U Decomposed Matrices
template <FPNumber T>
inline Matrix<T> LUSolve(const LTMatrix<T>& L, const UTMatrix<T>& U, const Matrix<T>& RHSMatrix) {Matrix<T> b(RHSMatrix); LUSolveInPlace(L,U,b); return b;}

// Solve by L-U Decomposed Matrix Compact form (forward + backward substitution) (in-place)
template <FPNumber T>
inline void LUSolveInPlace(const Matrix<T>& LU, Matrix<T>& RHSMatrix)
{
	std::uint64_t NumEqs = LU.Rows(); std::uint64_t NumRhs = RHSMatrix.Cols(); std::uint64_t NumUnKnowns = RHSMatrix.Rows(); std::int64_t NumEqsMinus1 = NumEqs-1;
	if (NumEqs != NumUnKnowns || NumEqs != LU.Rows() || NumEqs != LU.Cols()) {std::println("Dimensions mismatched... Quitting!"); return;}
	for (std::uint64_t r{1}; r<NumEqs; ++r)
		for (std::uint64_t c{}; c<r; ++c)
			for (std::uint64_t i{}; i<NumRhs; ++i)
				RHSMatrix[r,i] -= LU[r,c]*RHSMatrix[c,i];
	std::vector<T> UrrInv; UrrInv.reserve(NumEqs); for (std::uint64_t i{}; i<NumEqs; ++i) UrrInv.push_back(1.0/LU[i,i]);
	for (std::uint64_t i{}; i<NumRhs; ++i)
	{
		RHSMatrix[NumEqsMinus1,i] *= UrrInv[NumEqsMinus1];
	}
	for (std::int64_t r{NumEqsMinus1-1}; r>-1; --r)
	{
		for (std::int64_t c{r+1}; c<NumEqs; ++c)
		{
			for (std::uint64_t i{}; i<NumRhs; ++i)
			{
				RHSMatrix[r,i] -= LU[r,c]*RHSMatrix[c,i];
			}
		}
		for (std::uint64_t i{}; i<NumRhs; ++i)
		{
			RHSMatrix[r,i] *= UrrInv[r];
		}
	}
}

// Solve by L-U Decomposed Matrices
template <FPNumber T>
inline Matrix<T> LUSolve(const Matrix<T>& LU, const Matrix<T>& RHSMatrix) {Matrix<T> b(RHSMatrix); LUSolveInPlace(LU,b); return b;}

/*****************************
 *           ____      ____  *
 *          / / /     / / /  *
 *        /////     /////    *
 *      / / /     / / /      *
 *    /////     /////        *
 *  |||||     |||||          *
 *    \\\\\     \\\\\        *
 *      \ \ \     \ \ \      *
 *        \\\\\     \\\\\    *
 *          \_\_\     \_\_\  *
 *                           *
 *****************************/
template <Number T>
inline Matrix<T> Matrix<T>::Inverse() const
{
	auto LU = LUDecomposeCompact(*this);
	auto InverseThis(Matrix<T>::I(nRows));
	if (nRows!=nCols) {std::println("Cannot inverse via L-U decomposition method! Returning Identity"); return InverseThis;}
	else
	{
		LUSolveInPlace(LU,InverseThis);
		return InverseThis;
	}
}

template <Number T>
inline T Matrix<T>::Determinant() const
{
	auto LU = LUDecomposeCompact(*this);
	T result{};
	if (nRows!=nCols) {std::println("Cannot calculate... Quitting!"); return result;}
	else {for (std::uint64_t i{}; i<nRows; ++i) result += LU[i,i]; return result;}
}

/*********************************************************************************
 * _            _                 __________  _            _                     *
 * |\          /|       /\       |____  ____| |\          /| ||       || ||      *
 * |\\        //|      //\\           ||      |\\        //| ||       || ||      *
 * ||\\      //||     //  \\          ||      ||\\      //|| ||       || ||      *
 * || \\    // ||    //____\\         ||      || \\    // || ||       || ||      *
 * ||  \\  //  ||   /________\        ||      ||  \\  //  || ||       || ||      *
 * ||   \\//   ||  //        \\       ||      ||   \\//   || ||\_____/|| ||____  *
 * ||    \/    || //          \\      ||      ||    \/    ||  \_______/  |_____| *
 *                                                                               *
 *********************************************************************************/

template <Number T>
inline Matrix<T> MatMul(const LTMatrix<T>& L, const UTMatrix<T>& U, bool ForceMultiply=false)
{
	std::uint64_t mid{}; std::uint64_t nside = L.Rows(); std::uint64_t nside2 = U.Rows();
	if (!ForceMultiply)
	{
		Matrix<T> result(nside,nside);
		if (nside==nside2)
		{
			for (std::uint64_t r{}; r<nside; ++r)
			{
				for (std::uint64_t c{}; c<nside; ++c)
				{
					mid = std::min(r,c)+1;
					for (std::uint64_t m{}; m<mid; ++m)
					{
						result[r,c]+=L[r,m]*U[m,c];
					}
				}
			}
			return result;
		}
		else
		{
			std::println("Dimension mismatch {} != {}...\nReturning a matrix of zeros of shape ({},{})",nside,nside2,nside,nside);
			return result;
		}
	}
	else
	{
		nside = std::min(nside,nside2); Matrix<T> result(nside,nside);
		for (std::uint64_t r{}; r<nside; ++r)
		{
			for (std::uint64_t c{}; c<nside; ++c)
			{
				mid = std::min(r,c)+1;
				for (std::uint64_t m{}; m<mid; ++m)
				{
					result[r,c]+=L[r,m]*U[m,c];
				}
			}
		}
		return result;
	}
}

template <Number T>
inline Matrix<T> MatMul(const UTMatrix<T>& U, const LTMatrix<T>& L, bool ForceMultiply=false)
{
	std::uint64_t mid{}; std::uint64_t nside = L.Rows(); std::uint64_t nside2 = U.Rows();
	if (!ForceMultiply)
	{
		Matrix<T> result(nside,nside);
		if (nside==nside2)
		{
			for (std::uint64_t r{}; r<nside; ++r)
			{
				for (std::uint64_t c{}; c<nside; ++c)
				{
					mid = std::max(r,c);
					for (std::uint64_t m{mid}; m<nside; ++m)
					{
						result[r,c]+=U[r,m]*L[m,c];
					}
				}
			}
			return result;
		}
		else
		{
			std::println("Dimension mismatch {} != {}...\nReturning a matrix of zeros of shape ({},{})",nside,nside2,nside,nside);
			return result;
		}
	}
	else
	{
		nside = std::min(nside,nside2); Matrix<T> result(nside,nside);
		for (std::uint64_t r{}; r<nside; ++r)
		{
			for (std::uint64_t c{}; c<nside; ++c)
			{
				mid = std::max(r,c);
				for (std::uint64_t m{mid}; m<nside; ++m)
				{
					result[r,c]+=U[r,m]*L[m,c];
				}
			}
		}
		return result;
	}
}

template <Number T>
inline Matrix<T> MatMul(const Matrix<T>& m, const LTMatrix<T>& l, bool ForceMultiply=false)
{
	std::uint64_t row = m.Rows(); std::uint64_t col = l.Cols(); std::uint64_t mid = m.Cols();
	Matrix<T> result(row,col);
	if (!ForceMultiply)
	{
		if (mid==col)
		{
			for (std::uint64_t r{}; r<row; ++r)
			{
				for (std::uint64_t c{}; c<col; ++c)
				{
					for (std::uint64_t k{c}; k<mid; ++k)
					{
						result[r,c] += m[r,k]*l[k,c];
					}
				}
			}
			return result;
		}
		else
		{
			std::println("Dimension mismatch... Quitting!");
			return result;
		}
	}
	else
	{
		mid = std::min(mid,col);
		for (std::uint64_t r{}; r<row; ++r)
		{
			for (std::uint64_t c{}; c<col; ++c)
			{
				for (std::uint64_t k{c}; k<mid; ++k)
				{
					result[r,c] += m[r,k]*l[k,c];
				}
			}
		}
		return result;
	}
}

template <Number T>
inline Matrix<T> MatMul(const LTMatrix<T>& l, const Matrix<T>& m, bool ForceMultiply=false)
{
	std::uint64_t row = l.Rows(); std::uint64_t col = m.Cols(); std::uint64_t mid = m.Rows();
	Matrix<T> result(row,col);
	if (!ForceMultiply)
	{
		if (mid==row)
		{
			for (std::uint64_t r{}; r<row; ++r)
			{
				for (std::uint64_t c{}; c<col; ++c)
				{
					for (std::uint64_t k{}; k<=r; ++k)
					{
						result[r,c] += l[r,k]*m[k,c];
					}
				}
			}
			return result;
		}
		else
		{
			std::println("Dimension mismatch... Quitting!");
			return result;
		}
	}
	else
	{
		mid = std::min(mid,row);
		for (std::uint64_t r{}; r<row; ++r)
		{
			for (std::uint64_t c{}; c<col; ++c)
			{
				for (std::uint64_t k{}; k<=r; ++k)
				{
					result[r,c] += l[r,k]*m[k,c];
				}
			}
		}
		return result;
	}
}

template <Number T>
inline Matrix<T> MatMul(const Matrix<T>& m, const UTMatrix<T>& l, bool ForceMultiply=false)
{
	std::uint64_t row = m.Rows(); std::uint64_t col = l.Cols(); std::uint64_t mid = m.Cols();
	Matrix<T> result(row,col);
	if (!ForceMultiply)
	{
		if (mid==col)
		{
			for (std::uint64_t r{}; r<row; ++r)
			{
				for (std::uint64_t c{}; c<col; ++c)
				{
					for (std::uint64_t k{}; k<=c; ++k)
					{
						result[r,c] += m[r,k]*l[k,c];
					}
				}
			}
			return result;
		}
		else
		{
			std::println("Dimension mismatch... Quitting!");
			return result;
		}
	}
	else
	{
		mid = std::min(mid,col);
		for (std::uint64_t r{}; r<row; ++r)
		{
			for (std::uint64_t c{}; c<col; ++c)
			{
				for (std::uint64_t k{}; k<=c; ++k)
				{
					result[r,c] += m[r,k]*l[k,c];
				}
			}
		}
		return result;
	}
}

template <Number T>
inline Matrix<T> MatMul(const UTMatrix<T>& l, const Matrix<T>& m, bool ForceMultiply=false)
{
	std::uint64_t row = l.Rows(); std::uint64_t col = m.Cols(); std::uint64_t mid = m.Rows();
	Matrix<T> result(row,col);
	if (!ForceMultiply)
	{
		if (mid==row)
		{
			for (std::uint64_t r{}; r<row; ++r)
			{
				for (std::uint64_t c{}; c<col; ++c)
				{
					for (std::uint64_t k{r}; k<mid; ++k)
					{
						result[r,c] += l[r,k]*m[k,c];
					}
				}
			}
			return result;
		}
		else
		{
			std::println("Dimension mismatch... Quitting!");
			return result;
		}
	}
	else
	{
		mid = std::min(mid,row);
		for (std::uint64_t r{}; r<row; ++r)
		{
			for (std::uint64_t c{}; c<col; ++c)
			{
				for (std::uint64_t k{r}; k<mid; ++k)
				{
					result[r,c] += l[r,k]*m[k,c];
				}
			}
		}
		return result;
	}
}
