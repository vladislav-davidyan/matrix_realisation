#include <cmath>
#include <iostream>
#include "cl_matrix.h"

    Matrix::Matrix() 
    {   
        m = n = 0;
        A = nullptr;
    }

    Matrix::Matrix(int m_, int n_)
    {
        //std::cout<<"------------in Matrix initialization constructor------------"<<std::endl;
        m = m_;
        n = n_;
        if ((m < 1)||(n < 1))
        {
            std::cout << "Incorrect number of rows/columns"<< std::endl;
            A = nullptr;
        }
        else
        {
            A = (double**) new double*[m];

            for (int i = 0; i < m; ++i)
                A[i] = (double*) new double[n];

            for (int i = 0; i < m; ++i)
                for (int j = 0; j < n; ++j)
                    A[i][j] = 0.0;
        }
    }

    Matrix::Matrix(char* str)
    {
        //std::cout<<"------------in Matrix string initialization constructor------------"<<std::endl;
        int i=0;
        int num_of_strok=-1;
        while (str[i]!='\0')
        {    std::ostream& operator<<(std::ostream& c, Matrix& B);
        
            if (str[i]=='}') num_of_strok+=1;
            ++i;
        }
        int num_of_stolb=0;
        i=2;
        while (str[i]!='}')
        {
            if (str[i]==',') ++num_of_stolb;
            ++i;
        }
        if (num_of_stolb!=0) num_of_stolb+=1;
        //std::cout<<"row number "<<num_of_strok<<" "<<std::endl;
        //std::cout<<"column number "<<num_of_stolb<<" "<<std::endl;
        m=num_of_strok;
        n=num_of_stolb;
        i=0;
        int j=0;
        int k=0;
        int p=0;
        int st=0;
        int minus=0;
        double num=0;
        A=(double**) new double*[m];
        for (int i=0; i<m; i++)
            A[i] = (double*) new double[n];
        for (i=0; i<m; i++)
        {
            for (j=0; j<n; j++)
	        {
                while(str[k]!=',')
	            {
		            if ((i==m-1)&&(str[k]=='}')) break;
                    if ((str[k]>='0')&&(str[k]<='9')&&(p==0))
                        num=num*10+(str[k]-'0');
                    
	                if ((str[k]>='0')&&(str[k]<='9')&&(p==1))
                    {
                        num=num*10+(str[k]-'0');
                        ++st;
                    }
	                if (str[k]=='.') p=1;
		            if (str[k]=='-') minus=1;
		            ++k;
                }
	            A[i][j]=pow(-1,minus)*((double)(num/pow(10,st)));
                //std::cout<<A[i][j]<<std::endl;
                p=0;
	            st=0;
	            num=0;
	            minus=0;
	            ++k;
	        }
        }

    }

    Matrix::Matrix(double value)
    {
        std::cout<<"------------in 1-element Matrix constructor------------"<<std::endl;
        m = 1;
        n = 1;
        A = (double**) new double**;
        A[0] = (double*) new double*;
        A[0][0] = value;
    }

    Matrix::Matrix(const Matrix& A_)
    {
        //std::cout<<"------------in copy constructor------------"<<std::endl;
        m = A_.m;
        n = A_.n;

        A = (double**) new double*[m];

        for (int i = 0; i < m; ++i)
            A[i] = (double*) new double[n];
        
        for (int i = 0; i < m; ++i)
            for(int j = 0; j < m; ++j)
                A[i][j] = A_.A[i][j];
    }

    Matrix::Matrix(double* vector, int n_)
    {
        std::cout<<"------------in row Matrix constructor------------"<<std::endl;
        m = 1;
        n = n_;

        A = (double**) new double*[m];
        A[0] = (double*) new double[n];

        for(int i = 0; i < n; ++i)
            A[0][i] = vector[i];
    }

    Matrix::Matrix(int m_, double* vector)
    {
        std::cout<<"------------in column Matrix constructor------------"<<std::endl;
        m = m_;
        n = 1;

        A = (double**) new double*[m];

        for (int i = 0; i < m; ++i)
            A[i] = (double*) new double*;

        for(int i = 0; i < m; ++i)
            A[i][0] = vector[i];
    }

    double Matrix::get_elem(int i_, int j_)
    {
        std::cout<<"------------in get_elem------------"<<std::endl;
        if (((i_ >= 0) || (i_<m)) && ((j_ < n) && (j_ >=0)))
            return A[i_][j_];
        else
        {
            std::cout << "There is no way to reach element with ["<<i_<<";"<<j_<<"] coordinates" << std::endl;
            return 0.0;
        }
    }


    int Matrix::rows()
    {
        return m;
    }

    int Matrix::columns()
    {
        return n;
    }

    Matrix::~Matrix()
    { 
        //std::cout<<"------------in destructor------------"<<std::endl;
        if (n > 0)
        {
            for (int i = 0; i < m; ++i)
                delete[] A[i];
        }

        if (m > 0) delete[] A;
    }


    Matrix Matrix::operator=(const Matrix& B_)
    {
        //std::cout<<"------------in operator = ------------"<<std::endl;
        if ((m != B_.m) || (n != B_.n))
        {
            
            if (n > 0)
            {
              for (int i = 1; i <= m; i++)
                delete[] A[i];
            }

            if (m > 0)
                delete[] A;

            m = B_.m;
            n = B_.n;

            A = (double**) new double*[m];
            for (int i = 1; i <= m; ++i)
                A[i] = (double*) new double[n];
        }

        for (int i = 0; i < m; ++i)
            for (int j = 0; j < n; ++j)
                A[i][j] = B_.A[i][j];
                
        //std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
        return *this;
    }

    Matrix Matrix::operator [](int num_i)
    {
        std::cout<<"------------in operator [] ------------"<<std::endl;
        --num_i;
        int num_j = 0;

        if (num_i + 1 <= this->m)
        {     
    	    num_j = this->n;
    	    double nm[num_j];

          	for (int i = 0; i<this->m; i++)
                for (int j = 0; j<this->n; j++)
                    if (i == num_i)
                        nm[j]=this->A[i][j];

    	    return Matrix(nm, num_j);
        }

        else if (num_i + 1 <= this->n)
        {
            num_j = this->m;
    	    double nm[num_j];

            for (int i = 0; i < this->m; ++i)
                for (int j = 0; j < this->n; ++j)
                    if (j == num_j)
                        nm[i]=this->A[i][j];
    	    
            return Matrix(num_j,nm);
        }

        else 
        {
            std::cout << "Incorrect number ("<< num_i <<") of rows/columns"<< std::endl;
    	    return Matrix();
        }
    } 

    Matrix Matrix::operator *(double val)
    {
        std::cout<<"------------in operator * ------------"<<std::endl;
        if (this->A)
        {
	        Matrix M = Matrix(*this);
            for (int i = 0; i < this->m; ++i)
                for (int j = 0; j < this->n; ++j)
                    M.A[i][j] = val*(this->A[i][j]);
       	    return M;
        }
        else return Matrix();
    }

    void Matrix::operator *=(double val)
    {
        std::cout<<"------------in operator *= ------------"<<std::endl;
        if (this->A)
        {
            for (int i = 0; i < this->m; ++i)
                for (int j = 0; j < this->n; ++j)
                    this->A[i][j] = val*(this->A[i][j]);
        }
        else return;
    }

    Matrix Matrix::operator +(Matrix& M2)
    {
        Matrix M1 = Matrix (*this);
        if ((M1.m == M2.m)&&(M1.n==M2.n))
        {
            for(int i = 0; i < this->m; ++i)
                for(int j = 0; j < this->n; ++j)
                {
                    M1.A[i][j] = M2.A[i][j] + this->A[i][j];
                    //std::cout<< M1.A[i][j]<<std::endl;
                }
        }
        else
        {
            std::cout<<"Matrices have different size!"<<std::endl;
        }
        return M1;
    }

    void Matrix::operator +=(Matrix& M2)
    {
        if ((this->m == M2.m)&&(this->n == M2.n))
        {
            for(int i = 0; i < this->m; ++i)
                for(int j = 0; j < this->n; ++j)
                {
                    this->A[i][j] = M2.A[i][j] + this->A[i][j];
                    //std::cout<< M1.A[i][j]<<std::endl;
                }
        }
        else
            std::cout<<"Matrices have different size!"<<std::endl;
        return;
    }    

    Matrix Matrix::operator -(Matrix& M2)
    {
        Matrix M1 = Matrix (*this);
        if ((M1.m == M2.m)&&(M1.n==M2.n))
        {
            for(int i = 0; i < this->m; ++i)
                for(int j = 0; j < this->n; ++j)
                {
                    M1.A[i][j] = M2.A[i][j] - this->A[i][j];
                    //std::cout<< M1.A[i][j]<<std::endl;
                }
        }
        else
        {
            std::cout<<"Matrices have different size!"<<std::endl;
        }
        return M1;
    }

    void Matrix::operator -=(Matrix& M2)
    {
        if ((this->m == M2.m)&&(this->n == M2.n))
        {
            for(int i = 0; i < this->m; ++i)
                for(int j = 0; j < this->n; ++j)
                {
                    this->A[i][j] = M2.A[i][j] - this->A[i][j];
                    //std::cout<< M1.A[i][j]<<std::endl;
                }
        }
        else
            std::cout<<"Matrices have different size!"<<std::endl;
        return;
    }

    Matrix Matrix::operator *(Matrix& M2)
    {
        Matrix C = Matrix(*this);
        if (this->n == M2.m)
        {

            for(int i = 0; i < M2.m; i++)
                for(int j = 0; j < M2.n; j++)
                {
                    C.A[i][j] = 0;
                    for(int k = 0; k < M2.n; k++)
                        C.A[i][j] += (this->A[i][k]) * M2.A[k][j];
                }
        }
        else
        {
            std::cout<<"Matrices have unsuitable size!"<<std::endl;
            C.A = nullptr;
        }
        return C;
    }

    void Matrix::operator *=(Matrix& M2)
    {
        Matrix C = Matrix(*this);
        if (this->n == M2.m)
        {
            for(int i = 0; i < M2.m; i++)
                for(int j = 0; j < M2.n; j++)
                {
                    C.A[i][j] = 0;
                    for(int k = 0; k < M2.n; k++)
                    {
                        C.A[i][j] += (this->A[i][k]) * M2.A[k][j];
                    }
                    std::cout << C.A[i][j];
                    this->A[i][j] = C.A[i][j];
                }   
        }
        else
        {
            std::cout<<"Matrices have unsuitable size!"<<std::endl;
        }
        return ;
    }

    Matrix Matrix::operator -()
    {
        for(int i= 0; i < this->m; ++i)
            for(int j = 0; j<this->n; ++j)
                this->A[i][j] = - this->A[i][j];
        return *this;
    }

    int Matrix::operator ==(const Matrix& M2_)
    {
        if ((this->m == M2_.m) && (this->n == M2_.n))
        {
            for (int i=0; i < this->m; ++i)
                for (int j=0; j < this->n; ++j)
                {
                    if (abs(this->A[i][j] - M2_.A[i][j]) > EPS)
                        return 0;
                }
            return 1;
        }
        std::cout<<"Matrices are not equal - they have different size!"<<std::endl;
        return 0;
    }

    int Matrix::operator !=(const Matrix& M2_)
    {
        if (*this == M2_)
            return 0;
        else
            return 1;
    }

    Matrix Matrix:: operator |(const Matrix& M2_)
    {
        Matrix C(this->m, (this->n)+(M2_.n));
        if (this->m == M2_.m)
        {
            
            for (int i = 0; i < this->m; ++i)
            {
                for (int j = 0; j < this->n; ++j)
                    C.A[i][j] = this->A[i][j];
                for (int j = this->n; j < (this->n) +(M2_.n); ++j)
                    C.A[i][j] = M2_.A[i][j-(this->n)];
            }
            
            return C;
        }
        else
        {
            std::cout << "Unsuitable size"<< std::endl;
            C.A == nullptr;
        }
        return C;
    }

    Matrix Matrix::operator /(const Matrix& M2_)
    {
        Matrix C((this->m)+(M2_.m),this->n );
        if (this->n == M2_.n)
        {
            
            for (int i = 0; i < this->m; ++i)
            {
                for (int j = 0; j< this->n; ++j)
                    C.A[i][j] = this->A[i][j];
            }
            for (int i = this->m; i < (this->m) + (M2_.m); ++i)
            {
                for (int j = 0; j< this->n; ++j)
                    C.A[i][j] = this->A[i-(this->m)][j];
            }           
            return C;
        }
        else
        {
            std::cout << "Unsuitable size"<< std::endl;
            C.A == nullptr;
        }
        return C;
    }


    void Matrix::set_elem(int i_, int j_, double value)
    {
        //std::cout<<"------------in set_elem------------"<<std::endl;
        if (((i_ >= 0) || (i_<m)) && ((j_ < n) && (j_ >=0)))
            A[i_][j_] = value;
        else
        {
            std::cout << "There is no way to reach element with ["<<i_<<";"<<j_<<"] coordinates" << std::endl;
            return;
        }

    }

    std::ostream& operator<<(std::ostream& c, Matrix& B)
    {
        c << "~Output:" << std::endl;
        if (B.A == nullptr)
        {
	        c << "Empty Matrix" << std::endl;
	        return c;
        }
        for (int i = 0; i < B.m; i++)
        {
            for (int j = 0; j < B.n; j++)
                c << B.A[i][j] << " ";
	        c <<std::endl;
        }
        return c;
    }

    Matrix Matrix::identity(int n)
    {
        std::cout<<"------------in identity------------"<<std::endl;
        Matrix B;
        if (n < 1)
        {
            std::cout << "Incorrect number ("<< n <<") of rows/columns"<< std::endl;
            B.A = nullptr;
        }
        B.A = (double**) new double*[n];

        for (int i = 0; i < n; ++i)
            B.A[i] = (double*) new double[n];

        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
            {
                if (i != j)
                    B.A[i][j] = 0.0;
                else
                    B.A[i][j] = 1;
            }

        return B;
    }

    Matrix Matrix::diagonal(int n, double* vector)
    {
        std::cout<<"------------in diagonal------------"<<std::endl;
        Matrix B;
        if (n < 1)
        {
            std::cout << "Incorrect number ("<< n <<") of rows/columns"<< std::endl;
            B.A = nullptr;
        }
        else
        {
            B.n = B.m = n;
            B.A = (double**) new double*[n];

            for (int i = 0; i < n; ++i)
                B.A[i] = (double*) new double[n];

            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                {
                    if (i != j)
                        B.A[i][j] = 0.0;
                    else
                        B.A[i][j] = vector[i];
                }
        }
        return B;
    }
    Matrix Matrix::operator !()
    {
        Matrix M = *this;
        for (int i = 0; i < M.m; ++i)
            for (int j=0; j < M.n; ++j)
                M.A[i][j]=this->A[j][i];
        *this = M;
        std::cout<< M <<std::endl;
        double max[this->n];
        double temp;
        std::cout<<"Max numbers in columns:"<<std::endl;
        for (int i = 0; i < this->n; ++i)
        {
            temp = this->A[0][i];
            for (int j = 0; j < this->m; ++j)
            {
                if (this->A[j][i] > temp)
                    temp = this->A[j][i];
            }
            max[i] = temp;
            std::cout << max[i] << ' ';
        } 
        std::cout<<std::endl;
        return *this;
    }