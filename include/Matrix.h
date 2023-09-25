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
    Matrix(const Matrix& matrix) = default; //copy
    Matrix(const int column, const int row, const double number);
    Matrix identity(const int n); // create a n identity matrix
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
    Matrix operator-(Matrix B); 
    Matrix operator*(Matrix B); 
    Matrix operator+(Matrix B);
    Matrix operator*(const double& alpha);

    bool operator==(Matrix& B);
    const Matrix operator|(Matrix& B); 
    Matrix dot(Matrix A, Matrix B); 
    Matrix elementOperation(double (*func) (double ele)); 
    

    //linear symmetric 
    bool choleskyDecomp(Matrix& L); 
    bool gaussElimination(Matrix& A, Matrix& B, bool pivot_enabled);
    bool gaussElimination_with_LR_decomp(Matrix& A, Matrix& z); 
    bool qr(Matrix A, Matrix& q, Matrix& r); 
    bool orthogonalIteration(Matrix A, Matrix& eigen_vector, std::vector<double>& eigen_value); 
    const double det();
    


private:
    std::vector<std::vector<double>> value; 
    int _col;
    int _row;
};


namespace MatrixOperation
{
    Matrix forwardSubstitution(const Matrix& L, const Matrix& b); 
    Matrix backwardSubstitution(const Matrix& R, const Matrix& b); 
    Matrix operator*(const double& alpha, const Matrix& A);
};

