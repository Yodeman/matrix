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
![default_init](https://user-images.githubusercontent.com/59335237/152627580-54cb002a-488f-4a24-80d4-44fdce066a04.png)


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
![unit_matrix](https://user-images.githubusercontent.com/59335237/152627588-b8688d36-6f0a-415b-be44-382d05406fa2.png)


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
![random](https://user-images.githubusercontent.com/59335237/152627594-f37964e6-43fa-4cf2-83a4-9be03892d08a.png)


```c++
#include <iostream>
#include "../matrix/matrix.hpp"

int main()
{
	// some matrix operations
	Matrix<double> m = Matrix<double>::rand(10, 10);
	std::cout << "Matrix m:\n" << m << std::endl;
	std::cout << "m's determinant: " << m.determinant() << std::endl;
	std::cout << "m's inverse:\n" << m.inverse() << std::endl;
	std::cout << "m * m's inverse:\n" << m * m.inverse() << std::endl;

	return 0;
}
```

![operations](https://user-images.githubusercontent.com/59335237/152627615-e6870d30-cb8f-4a00-8c23-98506ec43529.png)

