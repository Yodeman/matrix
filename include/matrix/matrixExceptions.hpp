#pragma once

#include <stdexcept>

namespace my_matrix {
    struct MatrixZeroDivision: public std::runtime_error{
        MatrixZeroDivision(const std::string& err) : std::runtime_error{err}
        {}
    };

    struct MatrixInvalidIndexing: public std::runtime_error{
        MatrixInvalidIndexing(const std::string& err) : std::runtime_error{err}
        {}
    };

    struct MatrixInvalidDiagonal: public std::runtime_error{
        MatrixInvalidDiagonal(const std::string& err) : std::runtime_error{err}
        {}
    };

    struct MatrixEliminationError: public std::runtime_error{
        MatrixEliminationError(const std::string& err) : std::runtime_error{err}
        {}
    };

    struct MatrixDeterminantError: public std::runtime_error{
        MatrixDeterminantError(const std::string& err) : std::runtime_error{err}
        {}
    };

    struct MatrixInverseError: public std::runtime_error{
        MatrixInverseError(const std::string& err) : std::runtime_error{err}
        {}
    };
} // my_matrix
