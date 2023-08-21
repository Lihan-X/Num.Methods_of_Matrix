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

    Matrix::Matrix(Matrix& matrix)
    {
        value = matrix.getValue();
        _row = matrix.getRow(); 
        _col = matrix.getCol();
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

    const int Matrix::getRow()
    {
        return _row;
    }

    const int Matrix::getCol()
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

    const std::vector<std::vector<double>> Matrix::getValue()
    {
        return value; 
    }

    //basic operation

    std::vector<double>& Matrix::operator[](int n)
    {
        return value[n];
    }

    double& Matrix::operator()(const unsigned int row, const unsigned int col)
    {
        return value[row][col];
    }

    Matrix Matrix::transpose()
    {
        Matrix x = Matrix(value);
        for (int i = 0; i < getRow(); i++)
        {
            for (int j = 0; j < getCol(); j++)
            {
                x[j][i]=value[i][j];
            }
        }
        return x;
    }

    Matrix Matrix::operator+(Matrix B)
    {
        if ((getCol()==B.getCol()) && (getRow() == B.getRow()))
        {
            for (int i = 0; i < getRow(); i++)
            {
                for (int j = 0; j < getCol(); j++)
                {
                    B[i][j]=B[i][j]+value[i][j];
                }
            }
            return B;
        }
        else
            throw "A,B has different size";
    }

    Matrix Matrix::operator*(const double alpha)
    {
        Matrix result = *this;
        for (int i = 0; i < result.getRow(); i++)
        {
            for (int j = 0; j < result.getCol(); j++)
            {
                result[i][j]=result[i][j]*alpha;
            }
        }
        return result;
    }

    bool Matrix::operator==(Matrix& B)
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

    Matrix Matrix::dot(Matrix& A, Matrix& B)
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
                {
                    sum += A[i][a]*B[a][j];
                    
                }
                product[i][j] = sum; 
            }
        }
        return product;

    }


    //linear symmetric 
    bool Matrix::choleskyDecomp(Matrix& A, Matrix& L)
    {
        L = Matrix(A.getRow(), A.getCol(), 0); 
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
                s = A[i][j] - sum;
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

        return 1;
    }

    Matrix Matrix::vorwaertsEinsetzen(Matrix& L, Matrix& b)
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

    Matrix Matrix::rueckwaertsEinsetzen(Matrix& R, Matrix& b)
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

    bool Matrix::gaussElimination(Matrix& A, Matrix& B, bool pivot_enabled)
    {
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
            return false;

        return true;
    }

    bool Matrix::gaussElimination_with_LR_decomp(Matrix& A, Matrix& z)
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

    bool Matrix::QR_decomp(Matrix& A, Matrix& B)
    {
        int n = A.getRow();
        int q = B.getCol();
        double alpha, beta;
        Matrix D = Matrix(n, n, 0);
        for (int s = 0; s < n; s++)
        {
            double mu = 0;
            for (int i = s; i < n; i++)
                mu += A[i][s]*A[i][s];
            if (mu < Matrix::esp*Matrix::esp)
                return false;
            if (A[s][s] < 0)
                mu = sqrt(mu);
            else
                mu = -sqrt(mu);
            D[s][s] = mu;
            A[s][s] = A[s][s] - mu;
            alpha = -mu*A[s][s];

            for (int j = s+1; j < n; j++)
            {
                double sum = 0;
                for (int i = s; i < n; i++)
                    sum += A[i][s]*A[i][j];
                beta = sum / alpha;
                for(int i = 0; i < n; i++)
                    A[i][j] = A[i][j]-beta*A[i][s];
            }

            for (int k = 0; k < q; k++)
            {
                double sum = 0;
                for (int i = s; i < n; i++)
                    sum += A[i][s]*B[i][k];
                beta = sum / alpha;
                for (int i = s; i < n; i++)
                    B[i][k] = B[i][k]-beta*A[i][s];
            }

        }
        if (fabs(A[n-1][n-1]) <= Matrix::esp)
            return false;

        return true;

    }

    Matrix Matrix::inverseL(Matrix& L)
    {
        Matrix result = L;
        for (int i = 0; i < L.getRow(); i++)
        {
            result[i][i] = 1 / L[i][i];
            for (int j = 0; i <= i-1; j++)
            {
                double sum = 0;
                for (int k = j; k <= i-1; k++)
                {
                    sum += L[i][k]*result[k][j];
                }
                result[i][j] = -sum/L[i][i]; 
            }
        }
        return result;
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
}



#if TEST
//tests
using namespace MatrixAlgorithm;
int main()
{
    Matrix vec1 = Matrix({{1,0,0},{0,2,0},{0,0,3}});
    Matrix vec2 = Matrix({{1,2,4},{2,13,23},{4,23,77}});
    Matrix trans = vec1.transpose();
    Matrix y = y.dot(vec1, vec2);

    std::cout << y.toString() << std::endl; 
    


    
}
#endif
