#pragma hdrstop
#include <iostream>
#include "../include/Matrix.h"
#define TEST 1


namespace MatrixAlgorithm
{
    Matrix::Matrix(std::vector<std::vector<double>> value)
    {
        int last=value[0].size();
        for (auto it = value.begin(); it != value.end(); it++)
        {
            if (it->size() == last)
            {
                last = it->size();
                continue;
            }
            else
                throw "Matrix Init Exception";
        }
        this->value = value;
        _row = value.size();
        _col = last;

    }

    Matrix::Matrix(const int row, const int col,  const double number)
    {
        value.resize(row);
        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < col; j++)
            {
                value[i].push_back(number);
            }
        }
        this->_row = row;
        this->_col = col;

    }

    std::string Matrix::to_string()
    {
        std::string vector_text = "";
        vector_text += '[';
        for (int i = 0; i < value.size(); i++)
        {
            vector_text += '[';
            for (int j = 0; j < value[i].size(); j++)
            {
                vector_text += std::to_string(value[i][j]);
                if (j+1 != value[i].size())
                    vector_text += ',';
            }
            vector_text += ']';
            if (i+1 != value.size())
                    vector_text += ",\n";
        }
        vector_text += "]";
        return vector_text;
    }

    bool Matrix::is_matrix()
    {
        int last=value[0].size();
        for (auto it = value.begin(); it != value.end(); it++)
        {
            if (it->size() == last)
            {
                last = it->size();
                continue;
            }
            else
                return false;
        }
        return true;
    }

    int Matrix::row()
    {
        return _row;
    }

    int Matrix::col()
    {
        return _col;
    }

    bool Matrix::is_symmetric()
    {
        if (_col != _row)
            return false;
        Matrix trans = transpose(*this);
        if  ( trans == *this)
            return true;
        else
            return false;
    }

    //basic operation

    std::vector<double>& Matrix::operator[](int n)
    {
        return value[n];
    }

    Matrix transpose( Matrix& mat)
    {
        Matrix x = Matrix(mat.col(), mat.row(), 0);
        for (int i = 0; i < mat.row(); i++)
        {
            for (int j = 0; j < mat.col(); j++)
            {
                x[j][i]=mat[i][j];
            }
        }
        return x;
    }

    Matrix operator+(Matrix& value1, Matrix& value2)
    {
        Matrix product = value1;
        if ((value1.col()==value2.col()) && (value1.row() == value2.row()))
        {
            for (int i = 0; i < product.row(); i++)
            {
                for (int j = 0; j < product.col(); j++)
                {
                    product[i][j]=product[i][j]+value2[i][j];
                }
            }
            return product;
        }
        else
            throw "A,B has different size";
    }

    Matrix operator*(const double alpha, Matrix& A)
    {
        Matrix result = A;
        for (int i = 0; i < result.row(); i++)
        {
            for (int j = 0; j < result.col(); j++)
            {
                result[i][j]=result[i][j]*alpha;
            }
        }
        return result;
    }

    bool operator==(Matrix& A, Matrix& B)
    {
        for (int i = 0; i < A.row(); i++)
        {
            for (int j = 0; j < A.col(); j++)
            {
                if (A[i][j] == B[i][j])
                    continue;
                else
                    return false;
            }
        }
        return true;
    }

    //linear symmetric 
    Matrix cholesky_decomp(Matrix& A)
    {
        Matrix L = Matrix(A.row(), A.col(), 0); 
        double sum = 0;
        double s; 
        for (int i = 0; i < L.row(); i++)
        {
            for (int j = 0; j <= i; j++)
            {
                sum = 0;
                for (int k = 0; k < j; k++)
                {
                    sum += L[i][k]*L[j][k];
                }
                s = A[i][j] - sum;
                if (j == i)
                {
                    if (s < ESP)
                        throw "A nicht positiv definit";
                    L[i][i] = sqrt(s);
                }
                else
                {
                    L[i][j] = s/L[j][j];
                }
            }
        }

        return L;
    }

    Matrix vorwaerts_einsetzen(Matrix& L, Matrix& b)
    {
        Matrix x = Matrix(b.row(), 1, 0);
        for (int i = 0; i < b.row(); i++)
        {
            double sum = 0;
            for (int j = 0; j < i; j++)
            {
                sum += L[i][j]*x[j][0];
            }
            x[i][0] = (b[i][0] - sum)/L[i][i];
        }
        return x;
    }

    Matrix rueckwaerts_einsetzen(Matrix& R, Matrix& b)
    {
        int n = b.row();
        Matrix x = Matrix(b.row(), 1, 0);
        for (int i = n-1; i >= 0; i--)
        {
            double sum = 0;
            for (int j = i+1; j < n; j++)
            {
                sum += R[i][j]*x[j][0];
            }
            x[i][0] = (b[i][0] - sum)/R[i][i];
        }
        return x;
    }

    void gauss_elimination(Matrix& A, Matrix& B)
    {
        //pivot suche
        int n = A.row();
        int q = B.col();
        double p; //pivotzeilen
        double piv;
        double l;
        
        for (int s = 0; s < n-1; s++)
        {
            //pivot suche
            p = s;
            piv = fabs(A[s][s]);
            for (int i = s+1; i < n; i++)
            {
                if (fabs(A[i][s]) <= piv)
                    continue;
                else
                {
                    p = i;
                    piv = fabs(A[i][s]);
                }
            }
            //Zeilenvertauschung
            if (p == s)
                continue;
            else
            {
                for (int j = s; j < n; j++)
                {
                    l = A[s][j];
                    A[s][j] = A[p][j];
                    A[p][j] = l;
                }
                for (int k = 0; k < q; k++)
                {
                    l = B[s][k];
                    B[s][k] = B[p][k];
                    B[p][k] = l;
                }
            }
            //Elimination
            for (int i = s+1; i < n; i++)
            {
                l = A[i][s]/A[s][s];
                A[i][s] = 0;
                for (int j = s+1; j < n; j++)
                    A[i][j] = A[i][j] - l*A[s][j];
                for (int k = s+1; k < n; k++)
                    B[i][k] = B[i][k] - l*B[s][k];
            }

        }
    }
}



#if TEST
//tests

using namespace MatrixAlgorithm;
int main()
{
    Matrix vec1 = Matrix({{1,0,0},{0,2,0},{0,0,3}});
    Matrix vec2 = Matrix({{1,2,4},{2,13,23},{4,23,77}});
    Matrix trans = transpose(vec1);
    Matrix L = cholesky_decomp(vec2);
    Matrix A = Matrix({{1,1,1},{1,1.001,5},{1,2,2}});
    Matrix B = Matrix(3,1,1);
    B[1][0] = 2;
    std::cout << A.to_string() << std::endl;
    std::cout << B.to_string() << std::endl;
    
    gauss_elimination(A, B);
    std::cout << A.to_string() << std::endl;
    std::cout << B.to_string() << std::endl;
}
#endif
