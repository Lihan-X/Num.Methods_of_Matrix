#pragma once
#include <vector>
#include <string>
#include <math.h>
#define TWINCAT_LIB 0

#if TWINCAT_LIB == 1
#define fabs fabs_
#define cos cos_ 
#define sin sin_ 
#endif

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
    static constexpr double esp =1.11e-16;
    double getMaximumNorm(); 
    double getEuklischNorm(); 


    //basic operation
    double& operator()(const unsigned int row, const unsigned int col);
    Matrix transpose();
    Matrix operator+(Matrix B);
    Matrix operator-(Matrix B); 
    Matrix operator*(const double alpha);
    Matrix operator*(Matrix& B); 
    bool operator==(Matrix& B);
    const Matrix operator|(Matrix& B); 
    Matrix dot(Matrix A, Matrix B); 
    Matrix elementOperation(double (*func) (double ele)); 
    

    //linear symmetric 
    bool choleskyDecomp(Matrix& A, Matrix& L); 
    Matrix forwardSubstitution(Matrix& L, Matrix& b); 
    Matrix backwardSubstitution(Matrix& R, Matrix& b); 
    bool gaussElimination(Matrix& A, Matrix& B, bool pivot_enabled);
    bool gaussElimination_with_LR_decomp(Matrix& A, Matrix& z); 
    bool qr(Matrix A, Matrix& q, Matrix& r); 
    bool orthogonalIteration(Matrix A, Matrix& eigen_vector, std::vector<double>& eigen_value); 
    Matrix inverseL(Matrix& L); 
    const double det();
    


private:
    std::vector<std::vector<double>> value; 
    int _col;
    int _row;
};


