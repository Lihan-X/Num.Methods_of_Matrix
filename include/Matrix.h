#pragma once
#include <vector>
#include <string>
#include <math.h>
#include <cstring>
#define TWINCAT_LIB 0

#if TWINCAT_LIB == 1
#define fabs fabs_
#define cos cos_ 
#define sin sin_ 
#define sqrt sqrt_
#else
#include "../include/hresult.h"
#endif


namespace MatrixOperation
{

    
    class Matrix
    {
        //initialise matrix
    public:
        Matrix() {};
        Matrix(std::vector<std::vector<double>> value);
        Matrix(const Matrix& matrix) = default; //copy
        Matrix(const int column, const int row, const double number);
        Matrix identity(const int n) const; // create a n identity matrix
        std::string toString();
        bool isMatrix();
        const int getRow() const;
        const int getCol() const;
        bool isSymmetric(); 
        const std::vector<std::vector<double>> getValue() const; 
        static constexpr double esp =1.11e-16;
        double getMaximumNorm(); 
        double getEuklischNorm(); 


        //basic operation
        std::vector<double>& operator[](int n);
        const std::vector<double>& operator[](int n) const; 
        double& operator()(const unsigned int row, const unsigned int col);
        Matrix transpose();
        Matrix operator*(Matrix B); 
        Matrix operator*(const double& alpha);

        bool operator==(const Matrix& B) const;
        const Matrix operator|(Matrix& B); 
        Matrix dot(Matrix A, Matrix B); 
        Matrix elementOperation(double (*func) (double ele)); 
        

        //linear symmetric 
        HRESULT lu(Matrix& L, Matrix& R, Matrix& z) const;
        HRESULT choleskyDecomp(Matrix& L); 
        HRESULT gaussElimination(Matrix& A, Matrix& B, bool pivot_enabled);
        HRESULT qr(Matrix& q, Matrix& r); 
        const double det();
        


    private:
        std::vector<std::vector<double>> value; 
        int _col;
        int _row;
    };


    
    Matrix forwardSubstitution(const Matrix& L, const Matrix& b); 
    Matrix backwardSubstitution(const Matrix& R, const Matrix& b);
    Matrix LDUIteration(const Matrix& A, const Matrix& B);  
    Matrix operator*(const double& alpha, Matrix A);
    Matrix operator-(const Matrix& A, const Matrix& B); 
    Matrix operator+(const Matrix& A, const Matrix& B);
    Matrix operator+=(Matrix& A, const Matrix& B);
};

