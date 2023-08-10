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

    Matrix transpose(Matrix& mat)
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
    Matrix cholesky_zerlegung(Matrix& A)
    {
        Matrix L = Matrix(A.row(), A.col(), 0); 
        double sum = 0;
        int k = 0;
        for (int i = 0; i < L.row(); i++)
        {
            for (int j = 0; j <= i; j++)
            {
                if (j == i)
                {
                    for (k = 0; k < i; k++)
                        sum += L[i][k]*L[i][k];
                    L[i][i] = sqrt(A[i][i] - sum);
                    sum = 0;
                    k = 0;
                }
                else
                {
                    for (k = 0; k < j; k++)
                        sum += L[i][k]*L[j][k];
                    L[i][j] = (A[i][j]- sum)/L[j][j];
                    sum = 0;
                    k = 0;
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
}



#if TEST
//tests

using namespace MatrixAlgorithm;
int main()
{
    Matrix vec1 = Matrix({{1,0,0},{0,2,0},{0,0,3}});
    Matrix vec2 = Matrix({{1,2,4},{2,13,23},{4,23,77}});
    Matrix trans = transpose(vec1);
    //Matrix L = cholesky_zerlegung(vec2);
    Matrix R = Matrix({{1,2,2},{0,-1,-5},{0,0,-9}});
    Matrix b = Matrix(3,1,3);
    b[1][0] = -13;
    b[2][0] = -27;
    Matrix result = rueckwaerts_einsetzen(R,b);
    std::cout << result.to_string() << std::endl;
}
#endif
