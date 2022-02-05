#include <iostream>
#include "matrix/matrix.hpp"
#include "matrixConfig.h"

int main()
{
	std::cout << "\n\t\tmatrix " << matrix_VERSION_MAJOR << "."
			  << matrix_VERSION_MINOR << "." << matrix_VERSION_PATCH
			  << "\n\n";
	
    Matrix<int> m{{10, 20, 30}, {40, 50, 60}};
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
    
    std::cout << m << std::endl;
    std::cout << m2 << std::endl;
    std::cout << m*20 << std::endl;
	std::cout << m3 << std::endl;
	std::cout << m3/0.03 << std::endl;
	std::cout << m4 << std::endl;
	std::cout << m+2.4 << std::endl;
	std::cout << m/0.3 << std::endl;
	std::cout << m%2 << std::endl;
	std::cout << m*2.5 << std::endl;
	std::cout << m5 << std::endl;
	std::cout << m6 << std::endl;
	std::cout << m7 << std::endl;
	std::cout << m8 << std::endl;
	try{
		std::cout << m8/(1-1) << std::endl;
	}catch(std::exception& e){
		std::cerr << e.what() << std::endl;
	}
	std::cout << m9 << std::endl;
	std::cout << m10 << std::endl;

	std::cout << m9({std::slice(2, 8, 1), std::slice(2, 8, 1)}) << std::endl;
	std::cout << m10({std::slice(2, 8, 2), std::slice(2, 8, 2)}) << std::endl;

	m.transpose_();
	auto m11 = m7.transpose();
	auto m12 = m9.transpose();
	auto m14 = m13.transpose();
	std::cout << m << std::endl;
	std::cout << m11 << std::endl;
	std::cout << m12 << std::endl;
	std::cout << m9*m10 << std::endl;
	std::cout << m9+m10 << std::endl;
	std::cout << m9+m12 << std::endl;
	std::cout << m9-m12 << std::endl;
	std::cout << m7*m11 << std::endl;
	
	try{
		std::cout << m7 << std::endl;
		int d = m7.determinant();
		std::cout << d << std::endl;
		std::cout << m8 << std::endl;
		double d1 = m8.determinant();
		std::cout << d1 << std::endl;

	}
	catch(std::exception& e){
		std::cerr << e.what() << std::endl;
	}	
	//auto sm = Matrix<int>::rand(10, 10);
	auto sm = Matrix<int>(
			{
			    {48, 48, 97, 1, 21, 71, 54, 43, 38, 57},
				{31, 23, 39, 37, 15, 75, 10, 25, 31, 79},
				{8, 91, 99, 99, 11, 19, 45, 2, 21, 94},
				{96, 82, 30, 62, 6, 37, 84, 30, 52, 72},
				{94, 77, 6, 37, 68, 95, 72, 33, 13, 30},
				{8, 79, 79, 21, 42, 95, 98, 51, 67, 39},
				{31, 81, 56, 35, 25, 63, 43, 36, 59, 99},
				{56, 20, 41, 40, 69, 14, 40, 91, 78, 56},
				{100, 79, 70, 75, 7, 23, 51, 68, 6, 12},
				{86, 32, 61, 12, 29, 94, 28, 96, 10, 14}
		}
	);
	std::cout << "sm:\n" << sm << std::endl;
	std::cout << sm.determinant() << std::endl;
	std::cout << "sm inverse\n" << sm.inverse() << std::endl;
	std::cout << (sm.inverse()) * sm << std::endl;


    return 0;
}

