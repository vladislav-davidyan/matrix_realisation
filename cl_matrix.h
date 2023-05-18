#pragma once

class Matrix
{
private:
    double** A;     
    int m, n;       //rows and columns
    static const int EPS = 0.001;
public:
    Matrix();
    Matrix(int m_, int n_);
    Matrix(char* str);
    Matrix(double value);
    Matrix(const Matrix& A_);
    Matrix(double* vector, int n_);
    Matrix(int m_, double* vector);
    ~Matrix();

    double get_elem(int i_, int j_);
    void set_elem(int, int, double);
    friend std::ostream& operator<<(std::ostream&, Matrix&);
    int rows();
    int columns();
    static Matrix identity(int);
    static Matrix diagonal(int, double*);
    Matrix operator=(const Matrix& B_);
    Matrix operator [](int num_i);
    Matrix operator *(double val);
    void operator *=(double val);
    Matrix operator +(Matrix& M2);
    void operator +=(Matrix& M2);
    Matrix operator -(Matrix& M2);
    void operator -=(Matrix& M2);
    Matrix operator *(Matrix& M2);
    void operator *=(Matrix& M2);
    Matrix operator -();
    int operator ==(const Matrix& M2_);
    int operator !=(const Matrix& M2_);
    Matrix operator |(const Matrix& M2_);
    Matrix operator /(const Matrix& M2_);
    Matrix operator !();
    
};
    void set_elem(int i_, int j_, double value);
    Matrix identity(int n);
    Matrix diagonal(int n, double* vector);
    