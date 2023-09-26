#pragma hdrstop

#include "../include/Matrix.h"
#define TEST 1

namespace MatrixOperation
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
                return;
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

    Matrix Matrix::identity(const int n)
    {
        Matrix I = Matrix(n,n,0);
        for (int i = 0; i < n; i++)
            I[i][i] = 1;
        return I;
    }

    std::string Matrix::toString()
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

    bool Matrix::isMatrix()
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

    const int Matrix::getRow() const
    {
        return _row;
    }

    const int Matrix::getCol() const
    {
        return _col;
    }

    bool Matrix::isSymmetric()
    {
        if (_col != _row)
            return false;
        Matrix trans = transpose();
        if  ( trans == *this)
            return true;
        else
            return false;
    }

    const std::vector<std::vector<double>> Matrix::getValue() const
    {
        return value; 
    }

    double Matrix::getMaximumNorm()
    {
        double max = 0;
        for (auto it = value.begin(); it != value.end(); it++)
        {
            for(auto it1 = it->begin(); it1 != it->end(); it1++)
            {
                if (fabs(*it1) > max)
                {
                    max = *it1; 
                }

            }
        }

        return max; 
    }

    double Matrix::getEuklischNorm()
    {
        double norm = 0; 
        if (_col == 1)
        {
            for (int i = 0; i < _row; i++)
            {
                norm += pow(value[i][0], 2);
            }
        }
        return sqrt(norm); 
    }

    //basic operation

    std::vector<double>& Matrix::operator[](int n)
    {
        return value[n];
    }
    const std::vector<double>& Matrix::operator[](int n) const
    {
        return value[n];
    }

    double& Matrix::operator()(const unsigned int row, const unsigned int col)
    {
        return value[row][col];
    }

    Matrix Matrix::transpose()
    {
        Matrix x = Matrix(_col, _row, 0);
        for (int i = 0; i < _row; i++)
        {
            for (int j = 0; j < _col; j++)
            {
                x[j][i]=value[i][j];
            }
        }
        return x;
    }


    Matrix Matrix::operator*(Matrix B)
    {
        Matrix result; 
        result = dot(*this, B); 
        return result; 
    }


    Matrix Matrix::operator*(const double& alpha)
    {
        Matrix result = *this;
        for (auto it_row = value.begin(); it_row != value.end(); it_row++)
        {
            for (auto it = it_row->begin(); it != it_row->end(); it++)
                *it=(*it)*alpha;

        }
        return result;
    }

    bool Matrix::operator==(const Matrix& B) const
    {
        for (int i = 0; i < getRow(); i++)
        {
            for (int j = 0; j < getCol(); j++)
            {
                if (value[i][j] == B[i][j])
                    continue;
                else
                    return false;
            }
        }
        return true;
    }

    const Matrix Matrix::operator|(Matrix& B)
    {
        if (_row != B.getRow())
            throw; 
        else
        {
            for (int i = 0; i < _row; i++)
            {
                value[i].reserve(_col+B.getCol());
                for (int j = 0; j < B.getCol(); j++)
                    value[i].push_back(B[i][j]); 
            }
        }
        return *this; 
    }

    Matrix Matrix::dot(Matrix A, Matrix B)
    {
        int m = A.getRow();
        int n = B.getCol();
        if ((A.getCol()) != (B.getRow()))
            throw "not possible! "; 
        Matrix product = Matrix(m, n, 0);
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j<n; j++)
            {
                double sum = 0;
                for (int a = 0; a < A.getCol(); a++)
                    sum += A[i][a]*B[a][j];
                
                product[i][j] = sum; 
            }
        }
        return product;

    }

    Matrix Matrix::elementOperation(double (*func) (double ele))
    {
        Matrix result = *this; 
        for (auto it_row = value.begin(); it_row != value.end(); it_row++)
        {
            for (auto it = it_row->begin(); it != it_row->end(); it++)
                *it = func(*it);
        }
        return result; 
    }


    //linear symmetric 
    HRESULT Matrix::choleskyDecomp(Matrix& L)
    {
        HRESULT hr = S_OK;
        L = Matrix(getRow(), getCol(), 0); 
        double sum = 0;
        double s; 
        for (int i = 0; i < L.getRow(); i++)
        {
            for (int j = 0; j <= i; j++)
            {
                sum = 0;
                for (int k = 0; k < j; k++)
                {
                    sum += L[i][k]*L[j][k];
                }
                s = value[i][j] - sum;
                if (j == i)
                {
                    if (s < Matrix::esp)
                        return 0;
                    L[i][i] = sqrt(s);
                }
                else
                {
                    L[i][j] = s/L[j][j];
                }
            }
        }

        return hr;
    }

    HRESULT Matrix::gaussElimination(Matrix& A, Matrix& B, bool pivot_enabled)
    {
        HRESULT hr = S_OK;
        //pivot suche
        int n = A.getRow();
        int q = B.getCol();
        double p; //pivotzeilen
        double piv;
        double l;
        
        for (int s = 0; s < n-1; s++)
        {
            if (pivot_enabled)
            {
                //pivot suche
                p = s;
                piv = fabs(A[s][s]);
                for (int i = s+1; i < n; i++)
                {
                    if (fabs(A[i][s]) > piv)
                    {
                        p = i;
                        piv = fabs(A[i][s]);
                    }
                }
                if (piv <= Matrix::esp)
                    return 0;
                //Zeilenvertauschung
                if (p != s)
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
            }
            //Elimination
            for (int i = s+1; i < n; i++)
            {
                l = A[i][s]/A[s][s];
                A[i][s] = 0;
                for (int j = s+1; j < n; j++)
                    A[i][j] = A[i][j] - l*A[s][j];
                for (int k = 0; k < q; k++)
                    B[i][k] = B[i][k] - l*B[s][k];
            }

        }

        if (fabs(A[n-1][n-1]) <= Matrix::esp)
            hr = S_FALSE;

        return hr;
    }

    HRESULT Matrix::gaussElimination_with_LR_decomp(Matrix& A, Matrix& z)
    {
        int n = A.getRow();
        z = Matrix(n,1,0);
        for (int i = 0; i < n; i++)
            z[i][0] = i;
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
                if (fabs(A[i][s]) > piv)
                {
                    p = i;
                    piv = fabs(A[i][s]);
                }
            }
            if (piv <= Matrix::esp)
                return 0;
            //Zeilenvertauschung
            if (p != s)
            {
                int j;
                for (j = 0; j < n; j++)
                {
                    l = A[s][j];
                    A[s][j] = A[p][j];
                    A[p][j] = l;
                }
                    j = z[s][0];
                    z[s][0] = z[p][0]; 
                    z[p][0] = j;    
            }
            //Elimination
            for (int i = s+1; i < n; i++)
            {
                l = A[i][s]/A[s][s];
                A[i][s] = l;
                for (int j = s+1; j < n; j++)
                    A[i][j] = A[i][j] - l*A[s][j];
            }

        }

        if (fabs(A[n-1][n-1]) <= Matrix::esp)
            return false;

        return true;

    }


    HRESULT Matrix::qr(Matrix& q, Matrix& r)
    {
        
    }


    const double Matrix::det()
    {
        int n = _col; 
        double p; //pivotzeilen
        double piv;
        double l;
        
        for (int s = 0; s < n-1; s++)
        {
            //pivot suche
            p = s;
            piv = fabs(value[s][s]);
            for (int i = s+1; i < n; i++)
            {
                if (fabs(value[i][s]) > piv)
                {
                    p = i;
                    piv = fabs(value[i][s]);
                }
            }
            if (piv <= Matrix::esp)
                return 0;
            //Zeilenvertauschung
            if (p != s)
            {
                for (int j = s; j < n; j++)
                {
                    l = value[s][j];
                    value[s][j] = value[p][j];
                    value[p][j] = l;
                }
            }
            
            //Elimination
            for (int i = s+1; i < n; i++)
            {
                l = value[i][s]/value[s][s];
                value[i][s] = 0;
                for (int j = s+1; j < n; j++)
                    value[i][j] = value[i][j] - l*value[s][j];
            }

        }

        double d = 1;
        for (int i = 0; i < n; i++)
            d *= value[i][i]; 
        return d;
    }


    Matrix forwardSubstitution(const Matrix& L, const Matrix& b)
    {
        int q = b.getCol();
        int n = b.getRow();
        Matrix x = Matrix(n, q, 0);
        for (int k =0; k < q; k++)
        {
            for (int i = 0; i < n; i++)
            {
                double sum = 0;
                for (int j = 0; j <= i-1; j++)
                {
                    sum += L[i][j]*x[j][k];
                }
                x[i][k] = (b[i][k] - sum)/L[i][i];
            }
        }
        return x;
    } 
    Matrix backwardSubstitution(const Matrix& R, const Matrix& b)
    {
        int n = b.getRow();
        int q = b.getCol();
        Matrix x = Matrix(b.getRow(), b.getCol(), 0);
        for (int k =0; k < q; k++)
        {
            for (int i = n-1; i >= 0; i--)
            {
                double sum = 0;
                for (int j = i+1; j < n; j++)
                {
                    sum += R[i][j]*x[j][k];
                }
                x[i][k] = (b[i][k] - sum)/R[i][i];
            }
        }
        return x;
    }

    Matrix operator*(const double& alpha, Matrix A)
    {
        for (int i = 0; i < A.getRow(); i++)
        {
            for (int j = 0; j < A.getCol(); j++)
                A[i][j]=A[i][j]*alpha;

        }
        return A;
    }

    Matrix LDUIteration(const Matrix& A, const Matrix& B)
    {
        Matrix L, D_inverse, U, x, x_last;
        x = Matrix(A.getRow(), B.getCol(), 1);

        L = Matrix(A.getRow(), A.getCol(), 0); 
        D_inverse = Matrix(A.getRow(), A.getCol(), 0); 
        U = Matrix(A.getRow(), A.getCol(), 0); 
        for (int i = 0; i < A.getRow(); i++)
        {
            for (int j = 0; j < A.getCol(); j++)
            {
                if (i < j)
                    L[i][j] = A[i][j];
                if (i == j)
                    D_inverse[i][j] = 1/A[i][j]; 
                if (i > j)
                    U[i][j] = A[i][j]; 
            }
        }

        unsigned int n = 0; 
        while(n < 1000)
        {
            x_last = x; 
            x = D_inverse*(B - L*x_last - U*x);
            n++; 
        }
    }

    Matrix operator-(const Matrix& A, const Matrix& B)
    {
        Matrix result = Matrix(A);
        if ((A.getCol()==B.getCol()) && (A.getRow() == B.getRow()))
        {
            for (int i = 0; i < A.getRow(); i++)
            {
                for (int j = 0; j < A.getCol(); j++)
                {
                    result[i][j]=A[i][j]-B[i][j];
                }
            }
            return result;
        }
        else
            throw;
    }

    Matrix operator+(const Matrix& A, const Matrix& B)
    {
        Matrix result = Matrix(A);
        if ((A.getCol()==B.getCol()) && (A.getRow() == B.getRow()))
        {
            for (auto i = 0; i < A.getRow(); i++)
            {
                for (auto j = 0; j < A.getCol(); j++)
                {
                    result[i][j]=A[i][j]+B[i][j];
                }
            }
            return result;
        }
        else
            throw;
    }

    Matrix operator+=(const Matrix& A, const Matrix& B)
    {
        Matrix result = Matrix(A);
        if ((A.getCol()==B.getCol()) && (A.getRow() == B.getRow()))
        {
            for (auto i = 0; i < A.getRow(); i++)
            {
                for (auto j = 0; j < A.getCol(); j++)
                {
                    result[i][j]=A[i][j]+B[i][j];
                }
            }
            return result;
        }
        else
            throw;
    }
};

    



#if TEST
//tests
#include <iostream>
using namespace MatrixOperation;
int main()
{
    Matrix vec1 = Matrix({{1,0,0},{0,2,0},{0,0,3}});
    Matrix vec2 = Matrix({{0,3,1},{0,4,-2},{2,1,2}});
    Matrix A = Matrix({{4,-1,1},{9,-8,9},{11,-11,12}});
    Matrix result = vec1 * vec2;
    Matrix b = Matrix(3, 1, 1);
    std::cout << "b is " << b.toString() << std::endl; 
    Matrix x = MatrixOperation::backwardSubstitution(vec1,b);
    std::cout << x.toString() << std::endl;
    
    
}
#endif
