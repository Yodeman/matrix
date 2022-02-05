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
template<typename T>
Matrix<T>::Matrix(std::vector<T>& v)
	:rows{1}, cols{v.size()}, ndims{1}, elems{v}, stride{std::make_pair(0, 1)}
{}

template<typename T>
Matrix<T>::Matrix(const size_t r, const size_t c)
	:rows{r}, cols{c}, elems(r*c), stride{std::make_pair(cols, 1)}
{
    assert(rows > 0); //, "Invalid row dimension!!!");
    assert(cols > 0); //, "Invalid column dimension!!!");
    ndims = (rows>1 && cols > 1) ? 2 : 1;
}

template<typename T>
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

template<typename T>
	template<typename R>
		Matrix<T>::Matrix(const Matrix<R>& m)
			:ndims{m.ndim()}, rows{m.shape().first}, cols{m.shape().second}, offset{0}, stride{m.strides()}, elems(m.begin(), m.end()) {}

// element access
template<typename T>
T& Matrix<T>::operator()(const size_t& r, const size_t& c) const
{
    assert(r <= (rows-1) && c <= (cols-1));
    return elems[(offset + (stride.first*r) + (stride.second*c))];
    //return const_cast<T&>();
}

template<typename T>
Matrix<T> Matrix<T>::operator()(const std::array<std::slice, 2>& ind) const
{
	auto row_slice = ind[0];
	auto col_slice = ind[1];
	assert(row_slice.start()>= 0 && row_slice.size()>=0 && row_slice.stride()>0);
	assert(col_slice.start()>= 0 && col_slice.size()>=0 && col_slice.stride()>0);
	size_t rlength = static_cast<size_t>(ceil((row_slice.size() - row_slice.start()) / static_cast<double>(row_slice.stride())));
	size_t clength = static_cast<size_t>(ceil((col_slice.size() - col_slice.start()) / static_cast<double>(col_slice.stride())));
	if ((row_slice.size()==0) && (col_slice.size()==0)){
		return *this;
	}
	else if (clength==0) {
		Matrix<T> res(rlength, cols);
		for (auto i=0; i<rlength; ++i){
			auto r = row_slice.start() + (i * row_slice.stride());
			for (auto j=0; j<cols; ++j){
				res(i,j) = this->operator()(r,j);
			}
		}
		return std::move(res);
	}
	else if (rlength==0) {
		Matrix<T> res(rows, clength);
		for (auto i=0; i<rows; ++i){
			for (auto j=0; j<cols; ++j){
				auto c = col_slice.start() + (j * col_slice.stride());
				res(i,j) = this->operator()(i,c);
			}
		}
		return std::move(res);
	}
	else{
		Matrix<T> res(rlength, clength);
		for(auto i=0; i<rlength; ++i){
			auto r = row_slice.start() + (i * row_slice.stride());
			for (auto j=0; j<clength; ++j){
				auto c = col_slice.start() + (j * col_slice.stride());
				res(i,j) = this->operator()(r,c);
			}
		}
		return std::move(res);
	}
}

template<typename T>
std::vector<T> Matrix<T>::row(const size_t& i) const
{
	std::vector<T> v(cols);
	for (size_t j=0; j<cols; ++j){
	    v[j] = this->operator()(i, j);
	}
	return std::move(v);
}

template<typename T>
std::vector<T> Matrix<T>::col(const size_t& i) const
{
	std::vector<T> v(rows);
	for (size_t j=0; j<rows; ++j) {
	    v[j] = this->operator()(j, i);
	}
	return std::move(v);
}

// utilities
template<typename T>
    template<typename F>
		Matrix<T>& Matrix<T>::apply(F f)
{
	for(auto& x: elems)
	    f(x);
	return *this;
}

template<typename T>
    template<typename M, typename F>
		typename std::enable_if<std::is_same<Matrix<T>, M>::value, Matrix<T>&>::type Matrix<T>::apply(const M& m, F f)
{
	//for(auto i=begin(), j=m.begin(); i!=end(); ++i, ++j)
	for(size_t i=0; i<rows; ++i)
		for (size_t j=0; j<cols; ++j)
			f(this->operator()(i,j), m(i,j));
	return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator+=(const T& value)
{
	return apply([&](T& a){ a += value; });
}
template<typename T>
Matrix<T>& Matrix<T>::operator-=(const T& value)
{
	return apply([&](T& a){ a -= value; });
}
template<typename T>
Matrix<T>& Matrix<T>::operator*=(const T& value)
{
	return apply([&](T& a){ a *= value; });
}
template<typename T>
Matrix<T>& Matrix<T>::operator/=(const T& value)
{
	if(value==0) throw MatrixZeroDivision("operator/: zero division!!!");
	return apply([&](T& a){ a /= value; });
}
template<typename T>
Matrix<T>& Matrix<T>::operator%=(const T& value)
{
	return apply([&](T& a){ a %= value; });
}

template<typename T, typename R, typename RT=Common_type<T,R>>
Matrix<RT> operator+(Matrix<T>& m, const R& val)
{
	std::pair<size_t, size_t> m_shape = m.shape();
	Matrix<RT> res(m_shape.first, m_shape.second);
	for(size_t i=0; i<m_shape.first; ++i) {
		for(size_t j=0; j<m_shape.second; ++j) {
			res(i,j) = m(i,j) + val;
		}
	}
	return std::move(res);
}

template<typename T, typename R, typename RT=Common_type<T,R>>
Matrix<RT> operator-(Matrix<T>& m, const R& val)
{
	std::pair<size_t, size_t> m_shape = m.shape();
	Matrix<RT> res(m_shape.first, m_shape.second);
	for(size_t i=0; i<m_shape.first; ++i) {
		for(size_t j=0; j<m_shape.second; ++j) {
			res(i,j) = m(i,j) - val;
		}
	}
	return std::move(res);
}

template<typename T, typename R, typename RT=Common_type<T,R>>
Matrix<RT> operator*(Matrix<T>& m, const R& val)
{
	std::pair<size_t, size_t> m_shape = m.shape();
	Matrix<RT> res(m_shape.first, m_shape.second);
	for(size_t i=0; i<m_shape.first; ++i) {
		for(size_t j=0; j<m_shape.second; ++j) {
			res(i,j) = m(i,j) * val;
		}
	}
	return std::move(res);
}

template<typename T, typename R, typename RT=Common_type<T,R>>
Matrix<RT> operator/(Matrix<T>& m, const R& val)
{
	if(val==0) throw MatrixZeroDivision("operator/: zero division!!!");
	std::pair<size_t, size_t> m_shape = m.shape();
	Matrix<RT> res(m_shape.first, m_shape.second);
	for(size_t i=0; i<m_shape.first; ++i) {
		for(size_t j=0; j<m_shape.second; ++j) {
			res(i,j) = m(i,j) / val;
		}
	}
	return std::move(res);
}

template<typename T, typename R, typename RT=Common_type<T,R>>
Matrix<RT> operator%(Matrix<T>& m, const R& val)
{
	std::pair<size_t, size_t> m_shape = m.shape();
	Matrix<RT> res(m_shape.first, m_shape.second);
	for(size_t i=0; i<m_shape.first; ++i) {
		for(size_t j=0; j<m_shape.second; ++j) {
			res(i,j) = m(i,j) % val;
		}
	}
	return std::move(res);
}

template<typename T>
    template<typename M>
		typename std::enable_if<std::is_same<Matrix<T>,M>::value, Matrix<T>&>::type Matrix<T>::operator+=(const M& m)
{
	assert(m.rows==rows && m.cols==cols);
	return apply(m, [](T& a, const Value_type<M>& b){ a += b; });
}

template<typename T>
    template<typename M>
		typename std::enable_if<std::is_same<Matrix<T>, M>::value, Matrix<T>&>::type Matrix<T>::operator-=(const M& m)
{
	assert(m.rows==rows && m.cols==cols);
	return apply(m, [](T& a, const Value_type<M>& b){ a -= b; });
}

template<typename T, typename T2, typename RT = Common_type<T,T2>>
    Matrix<RT> operator+(const Matrix<T>& a, const Matrix<T2>& b)
{
	Matrix<RT> res = a;
	res += b;
	return std::move(res);
}

template<typename T, typename T2, typename RT = Common_type<T,T2>>
    Matrix<RT> operator-(const Matrix<T>& a, const Matrix<T2>& b)
{
	Matrix<RT> res = a;
	res -= b;
	return std::move(res);
}

template<typename T, typename T2>
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

template<typename T>
Matrix<T> Matrix<T>::rand(const size_t& r, const size_t& c)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	if(std::is_same<T,int>::value){
		std::uniform_int_distribution<> dis(1, 100);
		Matrix<T> res(r, c);
		for(auto& i : res)
			i = dis(gen);
		return std::move(res);
	}
	else if(std::is_same<T,double>::value){
		std::uniform_real_distribution<> dis(0, 10);
		Matrix<T> res(r, c);
		for(auto& i : res)
			i = dis(gen);
		return std::move(res);
	}
}

template<typename T>
Matrix<T> Matrix<T>::ones(const size_t& r, const size_t& c)
{
	Matrix<T> res(r, c);
	for(auto& i : res)
		i = 1;
	return std::move(res);
}

template<typename T>
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

template<typename T>
void Matrix<T>::transpose_()
{
	std::swap(stride.first, stride.second);
	std::swap(rows, cols);
}

template<typename T>
Matrix<T> Matrix<T>::transpose()
{
	Matrix<T> res(*this);
	res.transpose_();
	return std::move(res);
}

template<typename T>
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

template<typename T>
short Matrix<T>::elim_with_partial_pivot()
{
	short sign = 1;
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
		// wlimination
		for(size_t i=j+1; i<rows; ++i) {
			const double pivot = this->operator()(j,j);
			if(pivot == 0) throw MatrixEliminationError("can't solve: pivot==0");
			const double mult = this->operator()(i,j) / pivot;
			scale_and_add(j, i, mult);
		}
	}
	return sign;
}

template<typename T>
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

template<typename T>
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

template<typename T>
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

template<typename T>
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

template<typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& m)
{
	auto r = m.shape().first;
	auto c = m.shape().second;
	size_t rm=0, rn=0;
	T largest = m.elems.at(std::distance(m.cbegin(), std::max_element(m.cbegin(), m.cend(), [](const T& a, const T& b){return (std::abs(a) < std::abs(b));})));
	size_t width = std::to_string(largest).length() + 1;
	os.precision(4);
    os << "{\n";
    for (auto i=m.cbegin(); i!=m.cend(); ++i){
		if (rn==0) os << "    { ";
		//for(auto j : m.row(i)){
		os << std::right << std::setw(width) << ((std::abs(*i) > 1e-10) ? *i : 0) << " ";
		//}
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

