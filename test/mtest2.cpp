#include <iostream>
#include <matrix/matrix.hpp>
#include "matrixConfig.h"

int main()
{
	std::cout << "\n\t\tmatrix " << matrix_VERSION_MAJOR << "."
			  << matrix_VERSION_MINOR << "." << matrix_VERSION_PATCH
			  << "\n\n";
	
    Matrix<int> m{{10, 20, 30}, {40, 50, 60}};
    Matrix<int> m_0(2, 3);
    m_0 = m;
    Matrix<int> m_1(2, 3);
    Matrix<double> m2(4, 6);
    Matrix<double> m3{{3.14, 2.67, 9.1}, {3.8, 4.8, 6.5}};
    Matrix<int> m4 = m/3;
    auto m5 = Matrix<int>::unit(4, 4);
    auto m6 = Matrix<double>::unit(4, 4);
    auto m7 = Matrix<int>::rand(5, 5);
    auto m8 = Matrix<double>::rand(5, 5);
    auto m9 = Matrix<int>::rand(10, 10);
    auto m10 = Matrix<double>::rand(10,10);
    auto m13 = Matrix<int>::rand(1000, 3000);
    Matrix<int> m20(m3);
    Matrix<float> m21(m);
    Matrix<double> m23(m13);
    
    std::cout << "m : \n" << m << std::endl;
	std::cout << "m_0: \n" << m_0 << std::endl;
    std::cout << "m2 : \n" << m2 << std::endl;
    std::cout << "m*20 : \n" << m*20 << std::endl;
	std::cout << "m3 : \n" << m3 << std::endl;
	std::cout << "m3/0.03 : \n" << m3/0.03 << std::endl;
	std::cout << "m4 : \n" << m4 << std::endl;
	std::cout << "m+2.4 : \n" << m+2.4 << std::endl;
	std::cout << "m/0.3 : \n" << m/0.3 << std::endl;
	std::cout << "m%2 : \n" << m%2 << std::endl;
	std::cout << "m*2.5 : \n" << m*2.5 << std::endl;
	std::cout << "m5 : \n" << m5 << std::endl;
	std::cout << "m6 : \n" << m6 << std::endl;
	std::cout << "m7 : \n" << m7 << std::endl;
	std::cout << "m8 : \n" << m8 << std::endl;
	try{
		std::cout << "m8/(1-1) : \n" << m8/(1-1) << std::endl;
	}catch(std::exception& e){
		std::cerr << e.what() << std::endl;
	}
	std::cout << "m9 : \n" << m9 << std::endl;
	std::cout << "m10 : \n" << m10 << std::endl;

	std::cout << "m9[2:8:1, 2:8:1] : \n" << m9({std::slice(2, 8, 1), std::slice(2, 8, 1)}) << std::endl;
	std::cout << "m10[2:8:2, 2:8:2] : \n" << m10({std::slice(2, 8, 2), std::slice(2, 8, 2)}) << std::endl;
	std::cout << "m9[2:-1:1, 2:-1:1] : \n" << m9({std::slice(2, -1, 1), std::slice(2, -1, 1)}) << std::endl;
	std::cout << "m10[2:-1:2, 2:-1:1]" << m10({std::slice(2, -1, 2), std::slice(2, -1, 1)}) << std::endl;
	try{
		std::cout << "m9[-1:8:1, 1:8:1] : \n" << m9({std::slice(-1, 8, 1), std::slice(1, 8, 1)}) << std::endl;
		std::cout << "m9[2:8:1, 2:8:0]" << m9({std::slice(2, 8, 1), std::slice(2, 8, 0)}) << std::endl;
	}
	catch(std::exception& e){
		std::cerr << e.what() << std::endl;
	}

	m.transpose_();
	auto m11 = m7.transpose();
	auto m12 = m9.transpose();
	auto m14 = m13.transpose();
	std::cout << "transposed m : \n" << m << std::endl;
	std::cout << "m11 : \n" << m11 << std::endl;
	std::cout << "m12 : \n" << m12 << std::endl;
	std::cout << "m9*m10 : \n" << m9*m10 << std::endl;
	std::cout << "m9+m10 : \n" << m9+m10 << std::endl;
	std::cout << "m9+m12 : \n" << m9+m12 << std::endl;
	std::cout << "m9-m12 : \n" << m9-m12 << std::endl;
	std::cout << "m7*m11 : \n" << m7*m11 << std::endl;
	
	try{
		std::cout << "m7 : \n" << m7 << std::endl;
		int d = m7.determinant();
		std::cout << "m7.determinant() : \n" << d << std::endl;
		std::cout << "m8 : \n" << m8 << std::endl;
		double d1 = m8.determinant();
		std::cout << "m8.determinant() : \n" << d1 << std::endl;

	}
	catch(std::exception& e){
		std::cerr << e.what() << std::endl;
	}	
	auto sm = Matrix<int>::rand(10, 10);
	auto sm1 = Matrix<double>::rand(10, 10);
	std::cout << "sm:\n" << sm << std::endl;
	std::cout << sm.determinant() << std::endl;
	std::cout << "sm inverse\n" << sm.inverse() << std::endl;
	std::cout << (sm.inverse()) * sm << std::endl;

	std::cout << "sm1:\n" << sm1 << std::endl;
	std::cout << sm1.determinant() << std::endl;
	std::cout << "sm1 inverse\n" << sm1.inverse() << std::endl;
	std::cout << (sm1.inverse()) * sm1 << std::endl;

	m_1 = std::move(m_0);
	Matrix<double> m_2(std::move(m6));
	Matrix<double> m_3(5, 5);
	std::cout << "m_1 : \n" << m_1 << std::endl;
	std::cout << "m_2 : \n" << m_2 << std::endl;

	m_3 = m8;
	std::cout << "m_3 : \n" << m_3 << std::endl;
	std::cout << "m8 : \n" << m8 << std::endl;
	
	m_1 = m_1;
	std::cout << "self assignment\n";
	std::cout << "m_1 : \n" << m_1 << std::endl;

    return 0;
}

