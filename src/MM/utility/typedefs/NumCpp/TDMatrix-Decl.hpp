#pragma once
#include "Declarations.hpp"

// LTMatrix (Lower Triangular Matrix Class)
template <Number T>
class TDMatrix
{
	using value_type = T;
	private:
		std::uint64_t nDim{};
		std::vector<T> lo;
		std::vector<T> diag;
		std::vector<T> up;
		std::string name_of_the_type = get_type_name<value_type>();
		static constexpr const T zero{};
	public:
		/*******************
		 *                 *
		 *        X        *
		 *       /X\       *
		 *     ///X\\\     *
		 *   /////X\\\\\   *
		 *                 *
		 *******************/
		TDMatrix() = default;
		TDMatrix(const TDMatrix&) = default;
		TDMatrix(TDMatrix&&) noexcept = default;
		explicit TDMatrix(std::uint64_t n, T value=T{}) : nDim(n), lo((n>0)?n-1:0,value), diag(n,value), up((n>0)?n-1:0,value) {}
		TDMatrix(std::uint64_t n, T value_lo, T value_diag, T value_up) : nDim(n), lo((n>0)?n-1:0,value_lo), diag(n,value_diag), up((n>0)?n-1:0,value_up) {}
		TDMatrix(std::span<const T> span_lo, std::span<const T> span_diag, std::span<const T> span_up) : nDim(span_diag.size()), lo(span_lo.begin(),span_lo.end()), diag(span_diag.begin(),span_diag.end()), up(span_up.begin(),span_up.end()) {}
		/*********************************
		 *                               *
		 *   =========================   *
		 *   =========================   *
		 *                               *
		 *   =========================   *
		 *   =========================   *
		 *                               *
		 *********************************/
		TDMatrix& operator=(const TDMatrix&) = default;
		TDMatrix& operator=(TDMatrix&&) noexcept = default;
		TDMatrix& operator=(T value) {for (std::uint64_t i{}; i<nDim-1; ++i) {diag[i]=value; lo[i]=value; up[i]=value;} diag[nDim-1]=value; return *this;}
		/********************
		 *   ____    ____   *
		 *  |  __|  |__  |  *
		 *  | |        | |  *
		 *  | |        | |  *
		 *  | |        | |  *
		 *  | |        | |  *
		 *  | |__    __| |  *
		 *  |____|  |____|  *
		 *                  *
		 ********************/
		// operator[]
		// Matrix[i,j]
		inline T& operator[](std::uint64_t i, std::uint64_t j)
		{
			if (i>=nDim) throw std::runtime_error(std::to_string(i)+" out of bound!");
			if (j>=nDim) throw std::runtime_error(std::to_string(j)+" out of bound!");
			if (i==j) return diag[i];
			if (i==j+1) return lo[j];
			if (i==j-1) return up[i];
			throw std::runtime_error("Zero Element Write Access Denied");
		}
		inline const T& operator[](std::uint64_t i, std::uint64_t j) const
		{
			if (i>=nDim) throw std::runtime_error(std::to_string(i)+" out of bound!");
			if (j>=nDim) throw std::runtime_error(std::to_string(j)+" out of bound!");
			if (i==j) return diag[i];
			if (i==j+1) return lo[j];
			if (i==j-1) return up[i];
			return zero;
		}
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
		// info
		[[nodiscard]] inline std::uint64_t Rows() const noexcept { return nDim; }
		[[nodiscard]] inline std::uint64_t Cols() const noexcept { return nDim; }
		[[nodiscard]] inline std::uint64_t size() const noexcept { return 3*nDim-2; }
		[[nodiscard]] inline std::pair<std::uint64_t,std::uint64_t> Shape() const noexcept { return {nDim,nDim}; }
		[[nodiscard]] inline bool empty() const noexcept { return (lo.empty() && diag.empty() && up.empty()); }

		[[nodiscard]] inline std::vector<T>& GetLo() noexcept {return lo;}
		[[nodiscard]] inline std::vector<T>& GetDiag() noexcept {return diag;}
		[[nodiscard]] inline std::vector<T>& GetUp() noexcept {return up;}
		[[nodiscard]] inline const std::vector<T>& GetLo() const noexcept {return lo;}
		[[nodiscard]] inline const std::vector<T>& GetDiag() const noexcept {return diag;}
		[[nodiscard]] inline const std::vector<T>& GetUp() const noexcept {return up;}

		inline void printm(std::uint16_t width=10, std::uint16_t accuracy=8) const;
		inline Matrix<T> ToDense() const;
};

// print
template <Number T>
inline void TDMatrix<T>::printm(std::uint16_t width, std::uint16_t accuracy) const
{
	std::print("TDMatrix ({}), Shape: ({}, {})\n",name_of_the_type,nDim,nDim);
	if (nDim==0 || diag.empty()) return;
	std::print("{{");
	if constexpr (std::is_integral_v<T>)
	{
		for (size_t r=0; r<nDim; ++r)
		{
			r==0?std::print(" {{"):std::print("  {{");
			for (size_t c=0; c<nDim-1; ++c)
			{
				if (r==c) std::print(" {0:^{1}d},",diag[r],width);
				else if (r==c+1) std::print(" {0:^{1}d},",lo[c],width);
				else if (r==c-1) std::print(" {0:^{1}d},",up[r],width);
				else std::print(" {0:^{1}d},",T{},width);
			}
			if (r==nDim-1) std::print(" {0:^{1}d} }} }}\n",diag[r],width);
			else if (r==nDim-2) std::print(" {0:^{1}d} }}\n",up[r],width);
			else std::print(" {0:^{1}d} }}\n",T{},width);
		}
	}
	else
	{
		for (size_t r=0; r<nDim; ++r)
		{
			r==0?std::print(" {{"):std::print("  {{");
			for (size_t c=0; c<nDim-1; ++c)
			{
				if (r==c) std::print(" {0:^{1}.{2}g},",diag[r],width,accuracy);
				else if (r==c+1) std::print(" {0:^{1}.{2}g},",lo[c],width,accuracy);
				else if (r==c-1) std::print(" {0:^{1}.{2}g},",up[r],width,accuracy);
				else std::print(" {0:^{1}.{2}g},",T{},width,accuracy);
			}
			if (r==nDim-1) std::print(" {0:^{1}.{2}g} }} }}\n",diag[r],width,accuracy);
			else if (r==nDim-2) std::print(" {0:^{1}.{2}g} }}\n",up[r],width,accuracy);
			else std::print(" {0:^{1}.{2}g} }}\n",T{},width,accuracy);
		}
	}
}
