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
        Matrix() {};
        Matrix(std::vector<std::vector<double>> value);
        Matrix(Matrix& matrix); //copy
        Matrix(const int column, const int row, const double number);
        Matrix identity(const int n);
        std::string toString();
        bool isMatrix();
        const int getRow();
        const int getCol();
        std::vector<double>& operator[](int n);
        bool isSymmetric(); 
        const std::vector<std::vector<double>> getValue(); 
        static constexpr double esp =1e-10;
        double getMaximumNorm(); 
<<<<<<< HEAD
        double getEuklischNorm(); 
=======
        double getEuklidNorm(); 
>>>>>>> a9fb4ae56063c14c4e04d115c523f1be85c9ebd8


        //basic operation
        double& operator()(const unsigned int row, const unsigned int col);
        Matrix transpose();
        Matrix operator+(Matrix B);
        Matrix operator-(Matrix B); 
        Matrix operator*(const double alpha);
        bool operator==(Matrix& B);
        const Matrix operator|(Matrix& B); 
        Matrix dot(Matrix A, Matrix B); 
        

        //linear symmetric 
        bool choleskyDecomp(Matrix& A, Matrix& L); 
        Matrix vorwaertsEinsetzen(Matrix& L, Matrix& b); 
        Matrix rueckwaertsEinsetzen(Matrix& R, Matrix& b); 
        bool gaussElimination(Matrix& A, Matrix& B, bool pivot_enabled);
        bool gaussElimination_with_LR_decomp(Matrix& A, Matrix& z); 
        bool qr(Matrix A, Matrix& q, Matrix& r); 
        Matrix inverseL(Matrix& L); 
        const double det();
        


    private:
        std::vector<std::vector<double>> value; 
        int _col;
        int _row;
    };



}