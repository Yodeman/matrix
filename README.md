## matrix

### A simple implementation of 2D matrix.

###  (Work In Progress).

```c++
#include <iostream>
#include "matrix/matrix.hpp"

int main()
{
    // default initialization
    Matrix<int> m(10, 10);
    std::cout << m << std::endl;
    
    return 0;
}
```

```c++ 
#include <iostream>
#include "matrix/matrix.hpp"

int main()
{
    // unit matrix
    Matrix<int> m = Matrix<int>::unit(10, 10);
    std::cout << m << std::endl;
    
    return 0;
}
```

```c++
#include <iostream>
#include "matrix/matrix.hpp"

int main()
{
    // random matrix
    Matrix<double> m = Matrix<double>::rand(10, 10);
    std::cout << m << std::endl;
    
    return 0;
}
```

```c++
#include <iostream>
#include "../matrix/matrix.hpp"

int main()
{
	Matrix<double> m = Matrix<double>::rand(10, 10);
	std::cout << "Matrix m:\n" << m << std::endl;
	std::cout << "m's determinant: " << m.determinant() << std::endl;
	std::cout << "m's inverse:\n" << m.inverse() << std::endl;
	std::cout << "m * m's inverse:\n" << m * m.inverse() << std::endl;

	return 0;
}
```

