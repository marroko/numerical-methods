#include "opt/NR/numerical_recipes.c/nrutil.h"
#include "opt/NR/numerical_recipes.c/nrutil.c"
//#include "/opt/NR/numerical_recipes.c/gaussj.c"
//#include "/opt/NR/numerical_recipes.c/ludcmp.c"
//#include "/opt/NR/numerical_recipes.c/lubksb.c"

#include <iostream>
#include <string>
#include <iomanip>
#include <stdlib.h>
#include <math.h>
#include <chrono>

template <typename T>
class Vector
{
public:

    friend std::ostream& operator<<(std::ostream& output, const Vector<T> &vec)
    {
        output << "Vector " << vec.name << "\n";

        for(int i=1; i<=vec.size; ++i)
        {
            output << vec[i] << " ";
            output << "\n";
        }
        output << "\n";

        return output;
    }

    Vector(const int s = 1, std::string n = "") : size(s),
                                                  name(n) 
                                                  {
                                                      allocateVector();
                                                  }
    virtual ~Vector()
    {
        freeVector();
    }

    Vector & operator=(Vector &other) 
    {
        if(this == &other)
            return *this;

        if(this->size != other.getSize()) 
        {
            std::cout << "Zle przypisanie Vectora!\n";
            return *this;
        }

        freeVector();
        allocateVector();

        for(int i = 1 ; i <= size; i++)
            this->v[i] = other[i];

        return *this;
    }

    Vector operator-(const Vector &other) 
    {
        Vector tmp(other.size);

        for(int i = 1 ; i <= size; i++)
            tmp[i] = this->v[i] - other[i];

        return tmp;
    }

    Vector operator+(const Vector &other) 
    {
        Vector tmp(other.size);

        for(int i = 1 ; i <= size; i++)
            tmp[i] = this->f[i] + other[i];

        return tmp;
    }

    Vector operator*(float scalar) // mnozenie wektora przez skalar
    {
        Vector tmp(this->size);

        for(int i = 1 ; i <= size; i++)
            tmp[i] = this->v[i] * scalar;

        return tmp;
    }

    T operator*(const Vector &other) // skalarne mnozenie wektorow
    { 
        T sum = 0.0;

        for(int i = 1 ; i <= size; i++)
            sum += this->v[i] * other[i];

        return sum;
    }

    friend std::ostream& operator<<(std::ostream& output, Vector &vec) 
    {
        output << "Vector " << vec.name << "\n";

        for(int i=1; i<=vec.size; ++i)
            output << vec[i] << "\n";

        output << "\n";

        return output;
    }

    T normVector() // norma euklidesowa wektora jako pierwiastek z jego skalarnego mnozenia przez siebie
    {
        return sqrt(*this * *this);
    }

    T * operator()() { return v; }
    T & operator*() { return *v; }
    const T & operator[](unsigned i) const { return v[i]; }
    T & operator[](unsigned i) { return v[i]; }

    virtual void allocateVector() {}
    virtual void freeVector() {}

    int getSize() const { return size; }
    void setName(std::string n) { name = n; }

protected:
    T *v = nullptr;
    int size;
    std::string name;
};

class IVector : public Vector<int>
{
public:
    IVector(const int s = 1, std::string n = "") : Vector(s, n)
                                                   {}
    void allocateVector()
    {
        this->v = ivector(1, size);
    }

    void freeVector()
    {
        free_ivector(this->v, 1, size);
    }
};

class FVector : public Vector<float>
{
public:

    FVector(const int s = 1, std::string n = "") : Vector(s, n)
                                                   {}

    void allocateVector()
    {
        this->v = vector(1, size);
    }

    void freeVector()
    {
        free_vector(this->v, 1, size);
    }
};

class DVector : public Vector<double>
{
public:

    DVector(const int s = 1, std::string n = "") : Vector(s, n)
                                                   {}
    void allocateVector()
    {
        this->v = dvector(1, size);
    }

    void freeVector()
    {
        free_dvector(this->v, 1, size);
    }
};

template<typename T>
class Matrix 
{
public:

    friend std::ostream& operator<<(std::ostream& output, const Matrix<T> &m)
    {
        output << "\t\tMacierz " << m.name << "\n";

        for(int i=1; i<=m.rows; ++i)
        {
            for(int j=1; j<=m.columns; ++j)
                output << std::setw(15) << m(i,j) << " ";

            output << "\n";
        }
        output << "\n";

        return output;
    }

    Matrix(int r = 1, int c = 1, std::string n = "") :  rows(r),
                                                        columns(c),
                                                        name(n) 
                                                        {
                                                            allocateMatrix();
                                                        }
    virtual ~Matrix()
    {
        freeMatrix();
    }

    Matrix & operator=(const Matrix &other)
    {
        if(rows != other.rows || columns != other.columns) 
        {
            std::cout << "Zle przypisanie Matrixa!\n";
            return *this;
        }

        freeMatrix();
        allocateMatrix();

        for(int i = 1 ; i <= rows; i++) 
        {
            for(int j = 1 ; j <= columns; j++)
                this->m[i][j] = other.m[i][j];
        }
        return *this;
    }

    Matrix operator*(const Matrix& other) 
    {
        int n = this->rows;
        int m = other.columns;

        Matrix tmp(n,m);

        for(int i = 1 ; i <= n; i++) 
        {
            for(int j=1; j <= m; ++j) 
            {
                tmp.m[i][j] = 0.0;

                for(int k=1; k<=this->columns; ++k)
                    tmp.m[i][j] += this->m[i][k] * other.m[k][j];
            }
        }
        return tmp;
    }

    Matrix operator*(float scalar) 
    {
        int n = this->rows;
        int m = this->columns;

        Matrix tmp(n,m);

        for(int i = 1 ; i <= n; i++) 
        {
            for(int j=1; j <= m; ++j)
                tmp.m[i][j] = this->m[i][j] * scalar;
        }

        return tmp;
    }

    Matrix operator-(const Matrix &other) 
    {
        int n = this->rows;
        int m = this->columns;

        Matrix tmp(n,m);

        for(int i=1; i<=n; ++i) 
        {
            for(int j=1; j<=m; ++j)
                tmp.m[i][j] = this->m[i][j] - other(i,j);
        }

        return tmp;
    }

    Vector<T> operator*(const Vector<T> &other) // mnozenie Matrix przez FVector
    {
        int n = other.getSize();

        Vector<T> tmp(n);

        for(int i=1; i <= n; ++i)
            tmp[i] = 0.0;

        for(int i = 1; i <= this->rows; i++) {

            for(int j = 1; j <= this->columns; ++j)
                tmp[i] += this->m[i][j] * other[j];
        }

        return tmp;
    }

    T & operator()(int i, int j) { return m[i][j]; }
    const T & operator()(int i, int j) const { return m[i][j]; }
    T ** operator()() { return m; }

    T ptrMatrix() // wskaznik uwarunkowania Matrix jako najwieksza bezwzgledna wartosc z jej wnetrza
    {
        T max = fabs(m[1][1]);

        for(int i=1; i<=rows; ++i) 
        {
            for(int j=1; j<=columns; ++j)
            {
                if(max < fabs(m[i][j]))
                    max = fabs(m[i][j]);
            }
        }

        return max;
    }

    Matrix transpose() 
    {
        int n = this->rows;
        int m = this->columns;

        Matrix tmp(n,m);

        for(int i=1; i<=n; ++i) 
        {
            for(int j=1; j<=m; ++j)

                if(i == j)
                    tmp.m[i][j] = this->m[i][j];
                else
                    tmp.m[i][j] = this->m[j][i];
        }

        return tmp;
    }

    virtual void allocateMatrix() {}
    virtual void freeMatrix() {}

    void setName(std::string n) { name = n; }

protected:

    T **m;
    int rows;
    int columns;
    std::string name;
};

class FMatrix : public Matrix<float>
{
public:
	using Matrix<float>::operator(int i, int j);
	using Matrix<float>::operator();

    FMatrix(int r = 1, int c = 1, std::string n = "") : Matrix(r, c, n)
                                                        {}

    void freeMatrix()
    {
        free_matrix(m, 1, rows, 1, columns);
    }

    void allocateMatrix()
    {
        m = matrix(1, rows, 1, columns);
    }
};

class DMatrix : public Matrix<double>
{
public:

	using Matrix<double>::operator(int i, int j);
	using Matrix<double>::operator();


    DMatrix(int r = 1, int c = 1, std::string n = "") : Matrix(r, c, n)
                                                        {}

    void freeMatrix()
    {
        free_dmatrix(m, 1, rows, 1, columns);
    }

    void allocateMatrix()
    {
        m = dmatrix(1, rows, 1, columns);
    }
};

