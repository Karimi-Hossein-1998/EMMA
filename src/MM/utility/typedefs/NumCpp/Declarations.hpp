#pragma once
#include <cstddef>
#include <algorithm>
#include <initializer_list>
#include <iterator>
#include <print>
#include <ranges>
#include <vector>
#include <string>
#include <concepts>
#include <execution>
#include <stdfloat>
#include <cstdint>
#include <type_traits>
#include <stdexcept>
#include <span>
#include <cmath>
#include <complex>

static constexpr const size_t BLOCKSIZE = 32; // Side length of the square block of matrix (contiguous memory chunck of BLOCKSIZE*BLOCKSIZE size) to try and make operation cache-friendly!

// complex type
template <typename T>
struct is_complex : std::false_type {};
template <typename T>
struct is_complex<std::complex<T>> : std::true_type {};
template <typename T>
inline constexpr bool is_complex_v = is_complex<T>::value;

template <typename T>
concept Complex = is_complex_v<T>;
template <typename T>
concept FPComplex = is_complex_v<T> && std::floating_point<typename T::value_type>;
template <typename T>
concept Real = std::is_arithmetic_v<T>;
template <typename T>
concept FPReal = std::floating_point<T>;
template <typename T>
concept Number = Real<T> || Complex<T>;
template <typename T>
concept FPNumber = FPReal<T> || FPComplex<T>;

template <typename R>
concept Range2D = std::ranges::forward_range<R> && std::ranges::forward_range<std::ranges::range_reference_t<R>>;

inline std::uint64_t wrap_index(std::int64_t i,std::uint64_t nIndex)
{
	if (nIndex>0)
	{
		std::int64_t rem = i%static_cast<std::int64_t>(nIndex);
		return static_cast<std::uint64_t>(rem<0?rem+nIndex:rem);
	}
	else {return 0;}

}


#if __has_include(<stdfloat>)
#  include <stdfloat>
#endif

template <typename T>
std::string get_type_name()
{
	// 1. Guard against non-numeric types
	if constexpr (!Number<T>) return "unknown";
	// 2. Handle Complex types recursively using T::value_type
	else if constexpr (Complex<T>) return "complex<" + get_type_name<typename T::value_type>() + ">";
	// 3. Handle non-complex (Real scalar) numbers
	else
	{
		// Standard Floating-Point Types
		if constexpr (std::same_as<T, float>)                       return "float";
		else if constexpr (std::same_as<T, double>)                 return "double";
		else if constexpr (std::same_as<T, long double>)            return "long double";

		// C++23 Extended Floating-Point Types
		#if HAS_STDFLOAT
		else if constexpr (std::same_as<T, std::float16_t>)         return "std::float16_t";
		else if constexpr (std::same_as<T, std::bfloat16_t>)        return "std::bfloat16_t";
		else if constexpr (std::same_as<T, std::float32_t>)         return "std::float32_t";
		else if constexpr (std::same_as<T, std::float64_t>)         return "std::float64_t";
		else if constexpr (std::same_as<T, std::float128_t>)        return "std::float128_t";
		#endif

		#if defined(__GNUC__)
		else if constexpr (std::same_as<T,_Float32>)                return "std::float32";
		else if constexpr (std::same_as<T,_Float64>)                return "std::float64";
		else if constexpr (std::same_as<T,_Float128>)               return "std::float128";
		#endif

		// Fixed-Width Signed Integers
		else if constexpr (std::same_as<T, std::int8_t>)            return "std::int8_t";
		else if constexpr (std::same_as<T, std::int16_t>)           return "std::int16_t";
		else if constexpr (std::same_as<T, std::int32_t>)           return "std::int32_t";
		else if constexpr (std::same_as<T, std::int64_t>)           return "std::int64_t";

		// Fixed-Width Unsigned Integers
		else if constexpr (std::same_as<T, std::uint8_t>)           return "std::uint8_t";
		else if constexpr (std::same_as<T, std::uint16_t>)          return "std::uint16_t";
		else if constexpr (std::same_as<T, std::uint32_t>)          return "std::uint32_t";
		else if constexpr (std::same_as<T, std::uint64_t>)          return "std::uint64_t";

		// Fundamental Integer Types
		else if constexpr (std::same_as<T, char>)                   return "char";
		else if constexpr (std::same_as<T, signed char>)            return "signed char";
		else if constexpr (std::same_as<T, unsigned char>)          return "unsigned char";
		else if constexpr (std::same_as<T, short>)                  return "short";
		else if constexpr (std::same_as<T, unsigned short>)         return "unsigned short";
		else if constexpr (std::same_as<T, int>)                    return "int";
		else if constexpr (std::same_as<T, unsigned int>)           return "unsigned int";
		else if constexpr (std::same_as<T, long>)                   return "long";
		else if constexpr (std::same_as<T, unsigned long>)          return "unsigned long";
		else if constexpr (std::same_as<T, long long>)              return "long long";
		else if constexpr (std::same_as<T, unsigned long long>)     return "unsigned long long";

		// Compiler Extension Integers
		#if defined(__SIZEOF_INT128__)
		else if constexpr (std::same_as<T, __int128_t>)               return "int128_t";
		else if constexpr (std::same_as<T, unsigned __int128>)        return "unsigned int128";
		#endif

		else return "unknown";
	}
}

template <typename T>
class Ranges
{
private:
	T start=0;
	T stop=0;
	T step=1;
	std::uint64_t size=0;
public:
	Ranges() = default;
	Ranges(const Ranges&) = default;
	Ranges(Ranges&&) noexcept = default;
	explicit Ranges(T start_, T stop_, T step_) : start(start_), stop(stop_), step(step_)
	{
		if (step==0 || (step>0 && stop<=start) || (step<0 && stop>=start)) size = 0;
		else
		{
			std::float64_t diff = static_cast<std::float64_t>(stop-start); std::float64_t dStep = static_cast<std::float64_t>(step);
			std::float64_t divided = diff/dStep;
			std::uint64_t newSize = static_cast<std::uint64_t>(divided);
			std::float64_t eps=std::max<std::float64_t>(1e-14,divided*std::numeric_limits<std::float64_t>::epsilon()*10.0);
			size = (divided-static_cast<double>(newSize))>eps?newSize+1:newSize;
		}
	}
	[[nodiscard]] inline T Start() const noexcept {return start;}
	[[nodiscard]] inline T Stop() const noexcept {return stop;}
	[[nodiscard]] inline T Step() const noexcept {return step;}
	[[nodiscard]] inline std::uint64_t Size() const noexcept {return size;}
	[[nodiscard]] inline T operator[](std::int64_t index) const noexcept
	{
		if(index<size && size>=0) return start+index*step;
		else return start+wrap_index(index,size)*step;
	}
	Ranges<T>& operator=(const Ranges&) = default;
	Ranges<T>& operator=(Ranges&&) noexcept = default;
};

using Slice = Ranges<std::int64_t>;
