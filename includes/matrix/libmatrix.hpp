#ifndef MY_MATRIX
#define MY_MATRIX

#include <initializer_list>
#include <vector>
#include <valarray>
#include <array>
#include <type_traits>
#include "./matrixExceptions.hpp"

template<typename M>
using Value_type = typename M::value_type;

template<typename... T>
using Common_type = typename std::common_type<T...>::type;

template<typename T>
concept arithmetic_type = std::is_arithmetic_v<T>;

template<arithmetic_type T>
class Matrix;

template<arithmetic_type T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& m);

template<arithmetic_type T>
class Matrix{
	public:
		using value_type = T;
		using iterator = typename std::vector<T>::iterator;
		using const_iterator = typename std::vector<T>::const_iterator;

		// Constructors
		Matrix<T>() = delete;

		template<arithmetic_type R>
		Matrix<T>(const Matrix<R>&);

		Matrix<T>(Matrix<T>&& m)
			: ndims{m.ndims}, rows{m.rows}, cols{m.cols}, offset{m.offset},
			stride{m.stride}, elems(std::move(m.elems)) {}

		Matrix<T>& operator=(Matrix<T>&& m) {
			if (&m != this) {
				assert(m.rows == rows);
				assert(m.cols == cols);
				offset = m.offset;
				ndims = m.ndims;
				rows = m.rows;
				cols = m.cols;
				stride = m.stride;
				elems = std::move(m.elems);
			}
			return *this;
		}

		Matrix<T>(const Matrix<T>& m)
			: ndims{m.ndims}, rows{m.rows}, cols{m.cols}, offset{m.offset},
			stride{m.stride}, elems(m.cbegin(), m.cend()){}

		Matrix<T>& operator=(const Matrix<T>& m) {
			if (&m != this) {
				assert(m.rows == rows);
				assert(m.cols == cols);
				offset = m.offset;
				ndims = m.ndims;
				rows = m.rows;
				cols = m.cols;
				stride = m.stride;
				std::copy(m.cbegin(), m.cend(), elems.begin());
			}
			return *this;
		}

		explicit Matrix<T>(const size_t row, const size_t col);
		Matrix<T>(std::initializer_list<std::initializer_list<T>> v);
		
		~Matrix() = default;

		// Element access
		T& operator() (const size_t&, const size_t&);
		const T& operator()(const size_t&, const size_t&) const;
		Matrix<T> operator()(const std::array<std::slice, 2>&) const;
		Matrix<T> operator[](const size_t& i) const { return std::move(Matrix<T>(row(i))); }
		
		// utilities
		static Matrix<T> rand(const size_t& r, const size_t& c);
		static Matrix<T> unit(const size_t& r, const size_t& c);
		static Matrix<T> ones(const size_t& r, const size_t& c);

		iterator begin() { return elems.begin(); }
		iterator end() { return elems.end(); }
		const_iterator cbegin() const { return elems.cbegin(); }
		const_iterator cend() const { return elems.cend(); }

		double determinant();
		Matrix<double> inverse();
		Matrix<T> transpose();
		void transpose_();
		short elim_with_partial_pivot();

		template<typename F>
			Matrix& apply(F f);

		template<typename M, typename F>
			typename std::enable_if<std::is_same<Matrix<T>, M>::value, Matrix<T>&>::type apply(const M& m, F f);

		Matrix<T>& operator+=(const T& value);
		Matrix<T>& operator-=(const T& value);
		Matrix<T>& operator*=(const T& value);
		Matrix<T>& operator/=(const T& value);
		Matrix<T>& operator%=(const T& value);

		template<typename M>
			typename std::enable_if<std::is_same<Matrix<T>, M>::value, Matrix<T>&>::type operator+=(const M&);
		template<typename M>
			typename std::enable_if<std::is_same<Matrix<T>, M>::value, Matrix<T>&>::type operator-=(const M&);
		template<arithmetic_type T1, arithmetic_type T2>
			friend Matrix<Common_type<T1,T2>> operator*(const Matrix<T1>&, const Matrix<T2>&);

		size_t ndim() const { return ndims; }
		std::pair<size_t, size_t> shape() const { return std::make_pair(rows, cols); }
		std::pair<size_t, size_t> strides() const { return stride; }
		bool is_square() const { return rows==cols; }
		size_t size() const { return elems.size(); }
		T* data() { return elems.data(); }
		const T* data() const { return elems.data(); }

		friend std::ostream& operator<< <> (std::ostream& os, const Matrix<T>& m);
	
	protected:
		Matrix(std::vector<T>&);
		std::vector<T> row(const size_t& i) const;
		std::vector<T> col(const size_t& i) const;
		void swap_row(const size_t& first, const size_t& second);
		void scale_and_add(const size_t& p_row, const size_t& m_row, const double& multiplier);
		Matrix<T> minor(const size_t& r, const size_t& c);

	private:
		size_t ndims, rows, cols;
		size_t offset=0;
		std::pair<size_t, size_t> stride;
		std::vector<T> elems;
};

#endif // MY_MATRIX
