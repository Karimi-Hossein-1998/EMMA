#pragma once
#include "Declarations.hpp"

// LTMatrix (Lower Triangular Matrix Class)
template <Number T>
class LTMatrix
{
	using value_type = T;
	private:
		std::uint64_t nDim{};
		std::vector<T> data;
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
		LTMatrix() = default;
		LTMatrix(const LTMatrix&) = default;
		LTMatrix(LTMatrix&&) noexcept = default;
		explicit LTMatrix(std::uint64_t n, T value=T{}) : nDim(n), data(n*(n+1)/2,value) {}
		LTMatrix(std::span<const T> span_)
		{
			std::uint64_t ssize{span_.size()};
			data(span_.begin(),span_.end());
			std::uint64_t rMin{}, rMax{1}; while(rMin*rMax<(ssize<<1)) {++rMin; ++rMax;} nDim = rMax;
			if (((rMin*rMax)>>1)!=ssize) {data.reserve((rMin*rMax)>>1); data.insert(data.end(),((rMin*rMax)>>1)-ssize,T{});}
		}
		/*********************************
		 *                               *
		 *   =========================   *
		 *   =========================   *
		 *                               *
		 *   =========================   *
		 *   =========================   *
		 *                               *
		 *********************************/
		LTMatrix& operator=(const LTMatrix&) = default;
		LTMatrix& operator=(LTMatrix&&) noexcept = default;
		LTMatrix& operator=(T value) {for (std::uint64_t i{}; i<data.size(); ++i) data[i]=value; return *this;}
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
			if (i<j) throw std::runtime_error("Cannot assign values to zero upper triangles!");
			else return data[((i*(i+1))>>1)+j];
		}
		inline const T& operator[](std::uint64_t i, std::uint64_t j) const
		{
			if (i>=nDim) throw std::runtime_error(std::to_string(i)+" out of bound!");
			if (j>=nDim) throw std::runtime_error(std::to_string(j)+" out of bound!");
			if (i<j) return zero;// if (i<j) throw std::runtime_error("Cannot assign values to zero upper triangles!");
			else return data[((i*(i+1))>>1)+j];
		}
		// Matrix[i]
		inline std::span<T> operator[](std::uint64_t i)
		{
			if (i>=nDim) throw std::runtime_error(std::to_string(i)+" is out of bound");
			std::span<T> spanData {data}; return spanData.subspan((i*(i+1))>>1,i+1);
		}
		inline std::span<const T> operator[](std::uint64_t i) const
		{
			if (i>=nDim) throw std::runtime_error(std::to_string(i)+" is out of bound");
			std::span<T> spanData {data}; return spanData.subspan((i*(i+1))>>1,i+1);
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
		[[nodiscard]] inline std::uint64_t size() const noexcept { return data.size(); }
		[[nodiscard]] inline std::pair<std::uint64_t,std::uint64_t> Shape() const noexcept { return {nDim,nDim}; }
		[[nodiscard]] inline bool empty() const noexcept { return data.empty(); }
		[[nodiscard]] inline T* ptr() noexcept { return data.data(); }
		[[nodiscard]] inline const T* ptr() const noexcept { return data.data(); }

		inline void printm(std::uint16_t width=10, std::uint16_t accuracy=8) const;
		inline Matrix<T> ToDense() const;
};

// print
template <Number T>
inline void LTMatrix<T>::printm(std::uint16_t width, std::uint16_t accuracy) const
{
	std::print("LTMatrix ({}), Shape: ({}, {})\n",name_of_the_type,nDim,nDim);
	if (nDim==0 || data.empty()) return;
	std::print("{{");
	if constexpr (std::is_integral_v<T>)
	{
		for (size_t r=0; r<nDim; ++r)
		{
			r==0?std::print(" {{"):std::print("  {{");
			for (size_t c=0; c<nDim-1; ++c)
			{
				std::print(" {0:^{1}d},",r<c?T{}:data[((r*(r+1))>>1)+c],width);
			}
			(r==(nDim-1))?std::print(" {0:^{1}d} }} }}\n",data[data.size()-1],width):std::print(" {0:^{1}d}}},\n",T{},width);
		}
	}
	else
	{
		for (size_t r=0; r<nDim; ++r)
		{
			r==0?std::print(" {{"):std::print("  {{");
			for (size_t c=0; c<nDim-1; ++c)
			{
				std::print(" {0:^{1}.{2}g},",r<c?T{}:data[((r*(r+1))>>1)+c],width,accuracy);
			}
			(r==(nDim-1))?std::print(" {0:^{1}.{2}g} }} }}\n",data[data.size()-1],width,accuracy):std::print(" {0:^{1}.{2}g} }},\n",T{},width,accuracy);
		}
	}
}

template <Number T>
inline LTMatrix<T> MatMul(const LTMatrix<T>& l1, const LTMatrix<T>& l2, bool ForceMultiply=false)
{
	std::uint64_t ndim = l1.Rows(); std::uint64_t ndim2 = l2.Rows();
	if (ForceMultiply)
	{
		LTMatrix<T> result(ndim);
		if(ndim==ndim2)
		{
			for (std::uint64_t r{}; r<ndim; ++r)
			{
				for (std::uint64_t c{}; c<=r; ++c)
				{
					for (std::uint64_t k{c}; k<=r; ++k)
					{
						result[r,c] += l1[r,k]*l2[k,c];
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
		ndim = std::min(ndim, ndim2);
		LTMatrix<T> result(ndim);
		for (std::uint64_t r{}; r<ndim; ++r)
		{
			for (std::uint64_t c{}; c<=r; ++c)
			{
				for (std::uint64_t k{c}; k<=r; ++k)
				{
					result[r,c] += l1[r,k]*l2[k,c];
				}
			}
		}
		return result;
	}
}
