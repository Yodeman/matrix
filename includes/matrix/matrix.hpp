#include <iostream>
#include <cassert>
#include <numeric>
#include <iomanip>
#include <algorithm>
#include <random>
#include <cmath>
#include <string>
#include "libmatrix.hpp"

// constructors
template<arithmetic_type T>
Matrix<T>::Matrix(std::vector<T>& v)
	:rows{1}, cols{v.size()}, ndims{1}, elems{v}, stride{std::make_pair(0, 1)}
{}

template<arithmetic_type T>
Matrix<T>::Matrix(const size_t r, const size_t c)
	:rows{r}, cols{c}, elems(r*c), stride{std::make_pair(cols, 1)}
{
    assert(rows > 0); //, "Invalid row dimension!!!");
    assert(cols > 0); //, "Invalid column dimension!!!");
    ndims = (rows>1 && cols > 1) ? 2 : 1;
}

template<arithmetic_type T>
Matrix<T>::Matrix(std::initializer_list<std::initializer_list<T>> v)
{
    size_t sz = 0;
    rows = v.size();
    for(auto i = v.begin(); i!=v.end(); ++i){
		if (i==v.begin()) sz = i->size();
		assert(i->size() == sz);
		elems.insert(elems.end(), i->begin(), i->end());
    }
    cols = sz;
    stride = std::make_pair(cols, 1);
    ndims = ((rows>1) && (cols>1)) ? 2 : 1;
}

// construction of matrix from matrix of another type.
template<arithmetic_type T>
	template<arithmetic_type R>
		Matrix<T>::Matrix(const Matrix<R>& m)
			:ndims{m.ndim()}, rows{m.shape().first}, cols{m.shape().second}, offset{0}, stride{m.strides()}, elems(m.cbegin(), m.cend()) {}

// element access
template<arithmetic_type T>
T& Matrix<T>::operator()(const size_t& r, const size_t& c)
{
    if((0 > r || r >= rows) || (0 > c || c >= cols))
    	throw MatrixInvalidIndexing("Index out of range!!!");
    return elems[(offset + (stride.first*r) + (stride.second*c))];
}

template<arithmetic_type T>
const T& Matrix<T>::operator()(const size_t& r, const size_t& c) const
{
    if((0 > r || r >= rows) || (0 > c || c >= cols))
    	throw MatrixInvalidIndexing("Index out of range!!!");
    return elems[(offset + (stride.first*r) + (stride.second*c))];
}

template<arithmetic_type T>
Matrix<T> Matrix<T>::operator()(const std::array<std::slice, 2>& ind) const
{
	auto row_start = ind[0].start();
	auto row_stop = ind[0].size();
	auto row_stride = ind[0].stride();
	auto col_start = ind[1].start();
	auto col_stop = ind[1].size();
	auto col_stride = ind[1].stride();
	size_t rlength = 0, clength = 0;

	if (row_stop == -1)
		row_stop = rows;
	if (col_stop == -1)
		col_stop = cols;

	if (row_start < 0 || row_stop > rows || row_stride <= 0)
		throw MatrixInvalidIndexing("Invalid row slice index!!!");
	if (col_start < 0 || col_stop > cols || col_stride <= 0)
		throw MatrixInvalidIndexing("Invalid column slice index!!!");

	if (((row_stop - row_start) > 0) && ((col_stop - col_start) > 0)){
		rlength = static_cast<size_t>(ceil((row_stop - row_start) / static_cast<double>(row_stride)));
		clength = static_cast<size_t>(ceil((col_stop - col_start) / static_cast<double>(col_stride)));
	}

	if ((rlength==0) || (clength==0)){
		Matrix<T> res(rlength, clength);
		return std::move(res);
	}
	// returns the specifed rows and cols
	else{
		Matrix<T> res(rlength, clength);
		for(auto i=0; i<rlength; ++i){
			auto r = row_start + (i * row_stride);
			for (auto j=0; j<clength; ++j){
				auto c = col_start + (j * col_stride);
				res(i,j) = this->operator()(r,c);
			}
		}
		return std::move(res);
	}
}

template<arithmetic_type T>
std::vector<T> Matrix<T>::row(const size_t& i) const
{
	std::vector<T> v(cols);
	for (size_t j=0; j<cols; ++j){
	    v[j] = this->operator()(i, j);
	}
	return std::move(v);
}

template<arithmetic_type T>
std::vector<T> Matrix<T>::col(const size_t& i) const
{
	std::vector<T> v(rows);
	for (size_t j=0; j<rows; ++j) {
	    v[j] = this->operator()(j, i);
	}
	return std::move(v);
}

// utilities
template<arithmetic_type T>
    template<typename F>
		Matrix<T>& Matrix<T>::apply(F f)
{
	for(auto& x: elems)
	    f(x);
	return *this;
}

template<arithmetic_type T>
    template<typename M, typename F>
		typename std::enable_if<std::is_same<Matrix<T>, M>::value, Matrix<T>&>::type Matrix<T>::apply(const M& m, F f)
{
	auto i = begin();
	auto j = m.cbegin();
	while (j!=m.cend()){
		f(*i, *j);
		++j;
	}
	return *this;
}

template<arithmetic_type T>
Matrix<T>& Matrix<T>::operator+=(const T& value)
{
	return apply([&](T& a){ a += value; });
}
template<arithmetic_type T>
Matrix<T>& Matrix<T>::operator-=(const T& value)
{
	return apply([&](T& a){ a -= value; });
}
template<arithmetic_type T>
Matrix<T>& Matrix<T>::operator*=(const T& value)
{
	return apply([&](T& a){ a *= value; });
}
template<arithmetic_type T>
Matrix<T>& Matrix<T>::operator/=(const T& value)
{
	if(value==0) throw MatrixZeroDivision("operator/: zero division!!!");
	return apply([&](T& a){ a /= value; });
}
template<arithmetic_type T>
Matrix<T>& Matrix<T>::operator%=(const T& value)
{
	return apply([&](T& a){ a %= value; });
}

template<arithmetic_type T, arithmetic_type R, typename RT=Common_type<T,R>>
Matrix<RT> operator+(Matrix<T>& m, const R& val)
{
	Matrix<RT> res(m);
	return res.apply([&](RT& a){ a += val; });
}

template<arithmetic_type T, arithmetic_type R, typename RT=Common_type<T,R>>
Matrix<RT> operator-(Matrix<T>& m, const R& val)
{
	Matrix<RT> res(m);
	return res.apply([&](RT& a){ a -= val; });
}

template<arithmetic_type T, arithmetic_type R, typename RT=Common_type<T,R>>
Matrix<RT> operator*(Matrix<T>& m, const R& val)
{
	Matrix<RT> res(m);
	return res.apply([&](RT& a){ a *= val; });
}

template<arithmetic_type T, arithmetic_type R, typename RT=Common_type<T,R>>
Matrix<RT> operator/(Matrix<T>& m, const R& val)
{
	if(val==0) throw MatrixZeroDivision("operator/: zero division!!!");
	Matrix<RT> res(m);
	return res.apply([&](RT& a){ a /= val; });
}

template<arithmetic_type T, arithmetic_type R, typename RT=Common_type<T,R>>
Matrix<RT> operator%(Matrix<T>& m, const R& val)
{
	Matrix<RT> res(m);
	return res.apply([&](RT& a){ a %= val; });
}

template<arithmetic_type T>
    template<typename M>
		typename std::enable_if<std::is_same<Matrix<T>,M>::value, Matrix<T>&>::type Matrix<T>::operator+=(const M& m)
{
	assert(shape().first==m.shape().first && shape().second==m.shape().second);
	return apply(m, [](T& a, const Value_type<M>& b){ a += b; });
}

template<arithmetic_type T>
    template<typename M>
		typename std::enable_if<std::is_same<Matrix<T>, M>::value, Matrix<T>&>::type Matrix<T>::operator-=(const M& m)
{
	assert(shape().first==m.shape().first && shape().second==m.shape().second);
	return apply(m, [](T& a, const Value_type<M>& b){ a -= b; });
}

template<arithmetic_type T, arithmetic_type T2, typename RT = Common_type<T,T2>>
    Matrix<RT> operator+(const Matrix<T>& a, const Matrix<T2>& b)
{
	Matrix<RT> res = a;
	res += b;
	return std::move(res);
}

template<arithmetic_type T, arithmetic_type T2, typename RT = Common_type<T,T2>>
    Matrix<RT> operator-(const Matrix<T>& a, const Matrix<T2>& b)
{
	Matrix<RT> res = a;
	res -= b;
	return std::move(res);
}

template<arithmetic_type T, arithmetic_type T2>
Matrix<Common_type<T,T2>> operator*(const Matrix<T>& m1, const Matrix<T2>& m2)
{
	assert(m1.cols == m2.rows);
    size_t n = m1.rows;
    size_t m = m2.cols;
    Matrix<Common_type<T,T2>> res(n, m);
	for(int i=0; i<n; ++i) {
		auto r1 = m1.row(i);
		for(int j=0; j<m; ++j) {
			auto c2 = m2.col(j);
			res(i,j) = std::inner_product(r1.begin(), r1.end(), c2.begin(), static_cast<Common_type<T, T2>>(0));
		}
	}
	return std::move(res);
}

// generates random matrix. The supported types includes
// integers and doubles.
template<arithmetic_type T>
Matrix<T> Matrix<T>::rand(const size_t& r, const size_t& c)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	if(std::is_integral_v<T>){
		std::uniform_int_distribution<> dis(1, 100);
		Matrix<T> res(r, c);
		for(auto& i : res)
			i = dis(gen);
		return std::move(res);
	}
	else if(std::is_floating_point_v<T>){
		std::uniform_real_distribution<> dis(0, 10);
		Matrix<T> res(r, c);
		for(auto& i : res)
			i = dis(gen);
		return std::move(res);
	}
}

// generates a matrix whose elements are all 1
template<arithmetic_type T>
Matrix<T> Matrix<T>::ones(const size_t& r, const size_t& c)
{
	Matrix<T> res(r, c);
	for(auto& i : res)
		i = 1;
	return std::move(res);
}

template<arithmetic_type T>
Matrix<T> Matrix<T>::unit(const size_t& r, const size_t& c)
{
	if (r==c){
		Matrix<T> res(r, c);
		for(auto i=0; i<r; ++i){
			res(i,i) = 1;
		}
		return std::move(res);
	}
	throw MatrixInvalidDiagonal("row != col; a unit matrix is diagonal");
}

// transpose the matrix inplace.
template<arithmetic_type T>
void Matrix<T>::transpose_()
{
	std::swap(stride.first, stride.second);
	std::swap(rows, cols);
}

// returns a transposed matrix, with original matrix
// unchanged.
template<arithmetic_type T>
Matrix<T> Matrix<T>::transpose()
{
	Matrix<T> res(*this);
	res.transpose_();
	return std::move(res);
}

template<arithmetic_type T>
void Matrix<T>::swap_row(const size_t& first, const size_t& second)
{
	auto row_1 = row(first);
	auto row_2 = row(second);
	auto iter = begin();
	iter += (first*cols);
	for(auto j=row_2.begin(); j!=row_2.end(); ++j, ++iter)
		*iter = *j;
	iter = begin();
	iter += (second*cols);
	for(auto j=row_1.begin(); j!=row_1.end(); ++j, ++iter)
		*iter = *j;
}

// reduces the lower triangular part of the matrix to zeros
// using the guassian elimination with partial pivot method.
template<arithmetic_type T>
short Matrix<T>::elim_with_partial_pivot()
{
	short sign = 1; // used to keep track of the sign of determinant
					// as swapping a row affects the sign of the matrix's determinant.
	for(size_t j=0; j<rows; ++j) {
		size_t pivot_row = j;
		// look for a suitable pivot
		for(size_t k=j+1; k<rows; ++k) {
			if (std::abs(this->operator()(k, j)) > std::abs(this->operator()(pivot_row, j)))
				pivot_row = k;
		}
		// swap the rows if a better pivot is found
		if(pivot_row != j){
			swap_row(j, pivot_row);
			sign ^= -2;
		}
		// elimination
		for(size_t i=j+1; i<rows; ++i) {
			const double pivot = this->operator()(j,j);
			if(pivot == 0) throw MatrixEliminationError("can't solve: pivot==0");
			const double mult = this->operator()(i,j) / pivot;
			scale_and_add(j, i, mult);
		}
	}
	return sign;
}

template<arithmetic_type T>
void Matrix<T>::scale_and_add(const size_t& p_row, const size_t& m_row, const double& multiplier)
{
	auto pivot_row = row(p_row);
	auto other_row = row(m_row);
	auto iter = begin();
	iter += (m_row*cols);
	for(auto i=other_row.begin(), j=pivot_row.begin(); i!=other_row.end(); ++i, ++j, ++iter) {
		*iter = static_cast<T>(*i - ((*j) * multiplier));
	}
}

template<arithmetic_type T>
Matrix<T> Matrix<T>::minor(const size_t& r, const size_t& c)
{
	Matrix<T> minor_r_c(rows-1, cols-1);
	for(size_t i=0, m=0; i<rows; ++i) {
		if(i==r) continue;
		for (size_t j=0, n=0; j<cols; ++j) {
			if(j==c) continue;
			minor_r_c(m, n) = this->operator()(i,j);
			++n;
		}
		++m;
	}
	return std::move(minor_r_c);
}

template<arithmetic_type T>
Matrix<double> Matrix<T>::inverse()
{
	if(!is_square()) throw MatrixInverseError("Inverse of a non-square matrix!!!");
	auto d = this->determinant();
	if (d==0) throw MatrixInverseError("Inverse of a singular matrix is undefined!!!");
	if(rows==2 && cols==2) {
		Matrix<T> m{
			{this->operator()(1,1), (-1)*this->operator()(1,0)},
			{(-1)*this->operator()(0,1), this->operator()(0,0)}
		};
		return m/d;
	}
	Matrix<double>inv(rows, cols);
	size_t m=0, n=0;
	for(auto i=inv.begin(); i!=inv.end(); ++i){
		*i = ((minor(m,n).determinant()) * (pow(-1.0, static_cast<double>(m+n)))) / d;
		++n;
		if ((n == cols)) {
			++m;
			n = 0;
		}
	}
	inv.transpose_();
	return std::move(inv);
}

template<arithmetic_type T>
double Matrix<T>::determinant()
{
	if(!is_square()) throw MatrixDeterminantError("Determinant of non-square matrix!!!");
	double det = 1;
	Matrix<double> r(*this);
	auto sign = r.elim_with_partial_pivot();
	for(size_t i=0; i<rows; ++i)
		det *= r(i,i);
	return det*sign;
}

template<arithmetic_type T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& m)
{
	auto r = m.shape().first;
	auto c = m.shape().second;
	size_t rm=0, rn=0;
	T largest = m.elems.at(std::distance(m.cbegin(), std::max_element(m.cbegin(), m.cend(), [](const T& a, const T& b){return (std::to_string(a).length() < std::to_string(b).length());})));
	size_t width = std::to_string(largest).length() + 1;
	os.precision(4);
    os << "{\n";
    for (auto i=m.cbegin(); i!=m.cend(); ++i){
		if (rn==0) os << "    { ";
		os << std::right << std::setw(width) << ((std::abs(*i) > 1e-10) ? *i : 0) << " ";
		++rn;
		if ((rn==c)) {
			os << "}\n";
			++rm;
			rn = 0;
		}
    }
    os << "}\n" << '(' << m.shape().first << 'x' << m.shape().second << ") matrix\n";
    return os;
}

