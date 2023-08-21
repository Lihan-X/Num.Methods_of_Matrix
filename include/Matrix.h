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
        std::string toString();
        bool isMatrix();
        const int getRow();
        const int getCol();
        std::vector<double>& operator[](int n);
        bool isSymmetric(); 
        const std::vector<std::vector<double>> getValue(); 
        static constexpr double esp =1e-10;

        //basic operation
        double& operator()(const unsigned int row, const unsigned int col);
        Matrix transpose();
        Matrix operator+(Matrix B);
        Matrix operator*(const double alpha);
        bool operator==(Matrix& B);
        Matrix dot(Matrix& A, Matrix& B); 

        //linear symmetric 
        bool choleskyDecomp(Matrix& A, Matrix& L); 
        Matrix vorwaertsEinsetzen(Matrix& L, Matrix& b); 
        Matrix rueckwaertsEinsetzen(Matrix& R, Matrix& b); 
        bool gaussElimination(Matrix& A, Matrix& B, bool pivot_enabled);
        bool gaussElimination_with_LR_decomp(Matrix& A, Matrix& z); 
        bool QR_decomp(Matrix& A, Matrix& B); 
        Matrix inverseL(Matrix& L); 
        const double det();


    private:
        std::vector<std::vector<double>> value; 
        int _col;
        int _row;
    };



}