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

    Matrix Matrix::operator-(Matrix B)
    {
        if ((getCol()==B.getCol()) && (getRow() == B.getRow()))
        {
            for (int i = 0; i < getRow(); i++)
            {
                for (int j = 0; j < getCol(); j++)
                {
                    B[i][j]=value[i][j]-B[i][j];
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

    Matrix Matrix::operator*(Matrix& B)
    {
        Matrix result; 
        result = dot(*this, B); 
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

    const Matrix Matrix::operator|(Matrix& B)
    {
        if (_row != B.getRow())
            throw "not possible! "; 
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


    bool Matrix::qr(Matrix A, Matrix& q, Matrix& r)
    {
        int n = A.getRow();
        q = identity(n); 
        Matrix v, u, e;
        double mu;
        for (int s = 0; s < n-1; s++)
        {
            double sum = 0;
            for (int i = s; i < n; i++)
                sum += A(i,s)*A(i,s);
            if (A(s, s) < 0)
                mu = -sqrt(fabs(sum));
            else
                mu = sqrt(sum);
            v = Matrix(n-s, 1, 0); 
            e = v; 
            e(0,0) = mu; 
            for (int i = s; i < n; i++)
                v(i-s, 0) = A(i, s) - e(i-s, 0);
            u = v*(1/v.getEuklischNorm());
            Matrix q_temp = identity(n-s);
            Matrix q_ge = identity(n); 
            q_temp = q_temp - dot(u, u.transpose())*2;
            for (int i = s; i < n; i++)
            {
                for (int j = s; j < n; j++)
                    q_ge[i][j] = q_temp(i-s, j-s);
            } 
            q = dot(q_ge, q);
            A = dot(q_ge, A);
            
        }
        q = q.transpose(); 
        r = A;
        return true;
    }

    bool Matrix::orthogonalIteration(Matrix A, Matrix& eigen_vector, std::vector<double>& eigen_value)
    {
        const int iter = 200; 
        int count = 0;
        const int n = A.getCol(); 
        Matrix Y, Q, R; 
        eigen_value.resize(n); 
        double eigen_value_prev[n]; 
        for (int i = 0; i < n; i++) 
        {
            eigen_value[i] = 0;
            eigen_value_prev[i] = 1; 
        }
        Q = identity(n); 
        while (fabs(eigen_value[0]-eigen_value_prev[0]) > esp)
        {
            count++;
            if (count > iter)
                return false;
            for (int i = 0; i < n; i++)
                eigen_value_prev[i] = eigen_value[i];
            Y = A*Q;
            qr(Y, Q, R);
            for (int i = 0; i < n; i++)
                eigen_value[i] = R(i,i);
        }
        eigen_vector = Q;
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

    Matrix qrWithShiftAndDeflation(Matrix A)
    {
        
    }

    
}



#if TEST
//tests
using namespace MatrixAlgorithm;
int main()
{
    Matrix vec1 = Matrix({{1,0,0},{0,2,0},{0,0,3}});
    Matrix vec2 = Matrix({{0,3,1},{0,4,-2},{2,1,2}});
    Matrix A = Matrix({{4,-1,1},{9,-8,9},{11,-11,12}});
    Matrix trans = vec1.transpose();
    Matrix y,z,q,r;
    std::vector<double> eigen_value;
    vec1.orthogonalIteration(A, y, eigen_value); 
    std::cout << eigen_value[0] << ',' << eigen_value[1] << ',' << eigen_value[2] << std::endl;
    


    
}
#endif
