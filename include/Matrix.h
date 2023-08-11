#pragma once
#include <vector>
#include <string>
#include <math.h>
namespace MatrixAlgorithm
{
    #define ESP 1e-7
    //initialise matrix
    class Matrix
    {
    public:
        Matrix() {};
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
    bool cholesky_decomp(Matrix& A, Matrix& L); 
    Matrix vorwaerts_einsetzen(Matrix& L, Matrix& b); 
    Matrix rueckwaerts_einsetzen(Matrix& R, Matrix& b); 
    void gauss_elimination(Matrix& A, Matrix& B);

}