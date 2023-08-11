#pragma once
#include <vector>
#include <string>
#include <math.h>
namespace MatrixAlgorithm
{
    //initialise matrix
    class Matrix
    {
    public:
        Matrix(std::vector<std::vector<double>> value);
        Matrix(const int column, const int row, const double number);
        std::string to_string();
        bool is_matrix();
        int row();
        int col();
        std::vector<double>& operator[](int n);
        bool is_symmetric(); 

    private:
        std::vector<std::vector<double>> value; 
        int _col;
        int _row;
    };

    //basic operation
    Matrix transpose( Matrix& value);
    Matrix operator+(Matrix& value1, Matrix& value2);
    Matrix operator*(const double alpha, Matrix& A);
    bool operator==(Matrix& A, Matrix& B);

    //linear symmetric 
    Matrix cholesky_decomp(Matrix& A); 
    Matrix vorwaerts_einsetzen(Matrix& L, Matrix& b); 
    Matrix rueckwaerts_einsetzen(Matrix& R, Matrix& b); 
    Matrix gauss_elimination(Matrix& A, Matrix& b);

}