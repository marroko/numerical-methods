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

class IVector {

public:

    IVector(const long s = 1, std::string n = "") :
                                            size(s),
                                            name(n) { v = ivector(1,s);}

    ~IVector() { free_ivector(this->v, 1, size); }

    int * operator()() { return v; }
    int & operator*() { return *v; }
    const int & operator[](unsigned i) const { return v[i]; }
    int & operator[](unsigned i) { return v[i]; }

    IVector & operator=(IVector &other) {

        if(this == &other)
            return *this;
        if(this->size != other.getSize()) {

            std::cout << "Zle przypisanie IVectora!\n";
            return *this;
        }

        free_ivector(this->v, 1, size);
        this->v = ivector(1, size);

        for(long i = 1 ; i <= size; i++)
            this->v[i] = other[i];

       return *this;
    }
    IVector(const IVector &other) {

        this->size = other.getSize();
        this->v = ivector(1, size);

        for(long i = 1 ; i <= size; i++)
            this->v[i] = other[i];
    }

    friend std::ostream& operator<<(std::ostream& output, IVector &vec) {

        output << "IVector " << vec.name << "\n";

        for(long i=1; i<=vec.size; ++i) {

            output << vec[i] << " ";
            output << "\n";
        }
        output << "\n";

        return output;
    }

    long getSize() const { return size; }
    void setName(std::string n) { name = n; }

private:

    int *v = nullptr;
    long size;
    std::string name;
};

class FVector {

public:

    FVector(const long s = 1, std::string n = "") :
                                size(s),
                                name(n) { f = vector(1,s); }

    ~FVector() { free_vector(this->f, 1, size); }

    float * operator()() { return f; }
    float & operator*() { return *f; }
    const float & operator[](long i) const { return f[i]; }
    float & operator[](long i) { return f[i]; }

    FVector & operator=(const FVector &other) {

        if(this == &other)
            return *this;

        if(this->size != other.getSize()) {

            std::cout << "Zle przypisanie FVectora!\n";
            return *this;
        }

        free_vector(this->f, 1, size);
        this->f = vector(1, size);

        for(long i = 1 ; i <= size; i++)
            this->f[i] = other[i];

       return *this;
    }
    FVector(const FVector &other) {

        this->size = other.getSize();
        this->f = vector(1, size);

        for(long i = 1 ; i <= size; i++)
            this->f[i] = other[i];
    }

    FVector operator-(const FVector &other) {

        FVector tmp(other.size);
        for(long i = 1 ; i <= size; i++)
            tmp[i] = this->f[i] - other[i];

        return tmp;
    }

    FVector operator+(const FVector &other) {

        FVector tmp(other.size);

        for(long i = 1 ; i <= size; i++)
            tmp[i] = this->f[i] + other[i];

        return tmp;
    }

    FVector operator*(float scalar) { // mnozenie wektora przez skalar

        FVector tmp(this->size);

        for(long i = 1 ; i <= size; i++)
            tmp[i] = this->f[i] * scalar;

        return tmp;
    }

    float operator*(const FVector &other) {  // skalarne mnozenie wektorow

        float sum = 0.0;

        for(long i = 1 ; i <= size; i++)
            sum += this->f[i] * other[i];

        return sum;
    }

    friend std::ostream& operator<<(std::ostream& output, FVector &vec) {

        output << "FVector " << vec.name << "\n";

        for(long i=1; i<=vec.size; ++i)
            output << vec[i] << "\n";

        output << "\n";

        return output;
    }

    float normVector() { // norma euklidesowa wektora jako pierwiastek z jego skalarnego mnozenia przez siebie

        return sqrt(*this * *this);
    }

    long getSize() const { return size; }
    void setName(std::string n) { name = n; }

private:

    float *f = nullptr;
    long size;
    std::string name;
};

class DVector {

public:

    DVector(const long s = 1, std::string n = "") :  size(s),
                                                    name(n) { d = dvector(1,s); }

    ~DVector() { free_dvector(this->d, 1, size); }

    double * operator()() { return d; }
    double & operator*() { return *d; }
    const double & operator[](long i) const { return d[i]; }
    double & operator[](long i) { return d[i]; }

    DVector & operator=(const DVector &other) {

        if(this == &other)
            return *this;

        if(this->size != other.getSize()) {

            std::cout << "Zle przypisanie DVectora!\n";
            return *this;
        }

        free_dvector(this->d, 1, size);
        this->d = dvector(1, size);

        for(long i = 1 ; i <= size; i++)
            this->d[i] = other[i];

        return *this;
    }
    DVector(const DVector &other) {

        this->size = other.getSize();
        this->d = dvector(1, size);

        for(long i = 1 ; i <= size; i++)
            this->d[i] = other[i];
    }

    DVector operator-(const DVector &other) {

        DVector tmp(other.size);
        for(long i = 1 ; i <= size; i++)
            tmp[i] = this->d[i] - other[i];

        return tmp;
    }

    DVector operator+(const DVector &other) {

        DVector tmp(other.size);

        for(long i = 1 ; i <= size; i++)
            tmp[i] = this->d[i] + other[i];

        return tmp;
    }

    DVector operator*(double scalar) { // mnozenie wektora double przez skalar

        DVector tmp(this->size);

        for(long i = 1 ; i <= size; i++)
            tmp[i] = this->d[i] * scalar;

        return tmp;
    }

    double operator*(const DVector &other) {  // skalarne mnozenie wektorow double

        double sum = 0.0;

        for(long i = 1 ; i <= size; i++)
            sum += this->d[i] * other[i];

        return sum;
    }

    friend std::ostream& operator<<(std::ostream& output, DVector &vec) {

        output << "FVector " << vec.name << "\n";

        for(long i=1; i<=vec.size; ++i)
            output << vec[i] << "\n";

        output << "\n";
        return output;
    }

    double normVector() { // norma euklidesowa wektora double jako pierwiastek z jego skalarnego mnozenia przez siebie

        return sqrt(*this * *this);
    }

    long getSize() const { return size; }
    void setName(std::string n) { name = n; }

private:

    double *d = nullptr;
    long size;
    std::string name;
};

class Matrix {

public:

    ~Matrix() { free_matrix(this->m, 1, rows, 1, columns); }
    Matrix(long r = 1, long c = 1, std::string n = "") :
                                                m { matrix(1, r, 1, c) },
                                                rows(r),
                                                columns(c),
                                                name(n) {}

    Matrix & operator=(const Matrix &other){

       if(rows != other.rows || columns != other.columns) {

           std::cout << "Zle przypisanie Matrixa!\n";
           return *this;
       }
       free_matrix(m, 1, rows, 1, columns);
       m = matrix(1, rows, 1, columns);

        for(long i = 1 ; i <= rows; i++) {
            for(long j = 1 ; j <= columns; j++)
                this->m[i][j] = other.m[i][j];
        }
       return *this;
    }

    Matrix(const Matrix &other){

        rows = other.rows;
        columns = other.columns;
        m = matrix(1, rows, 1, columns);

       for(long i = 1; i <= rows; i++) {
           for(long j = 1 ; j <= columns; j++)
               this->m[i][j] = other.m[i][j];
       }
    }


    friend std::ostream& operator<<(std::ostream& output, const Matrix &m) {

        output << "\t\tMacierz " << m.name << "\n";

        for(long i=1; i<=m.rows; ++i) {

            for(long j=1; j<=m.columns; ++j)
                output << std::setw(15) << m(i,j) << " ";
            output << "\n";
        }
        output << "\n";

        return output;
    }

    Matrix operator*(const Matrix& other) {

        long n = this->rows;
        long m = other.columns;

        Matrix tmp(n,m);

        for(long i = 1 ; i <= n; i++) {

            for(long j=1; j <= m; ++j) {

                tmp.m[i][j] = 0.0;
                for(long k=1; k<=this->columns; ++k)
                    tmp.m[i][j] += this->m[i][k] * other.m[k][j];
            }
        }
        return tmp;
    }
    Matrix operator*(float scalar) {

        long n = this->rows;
        long m = this->columns;

        Matrix tmp(n,m);

        for(long i = 1 ; i <= n; i++) {

            for(long j=1; j <= m; ++j)
                tmp.m[i][j] = this->m[i][j] * scalar;
        }
        return tmp;
    }

    Matrix operator-(const Matrix &other) {

        long n = this->rows;
        long m = this->columns;

        Matrix tmp(n,m);

        for(long i=1; i<=n; ++i) {

            for(long j=1; j<=m; ++j)
                tmp.m[i][j] = this->m[i][j] - other(i,j);
        }
        return tmp;
    }

    FVector operator*(const FVector &other) { // mnozenie Matrix przez FVector

        long n = other.getSize();

        FVector tmp(n);

        for(long i=1; i <= n; ++i)
            tmp[i] = 0.0;

        for(long i = 1; i <= this->rows; i++) {

            for(long j = 1; j <= this->columns; ++j)
                tmp[i] += this->m[i][j] * other[j];
        }
        return tmp;
    }

    float & operator()(long i, long j) { return m[i][j]; }
    const float & operator()(long i, long j) const { return m[i][j]; }
    float ** operator()() { return m; }

    float ptrMatrix() { // wskaznik uwarunkowania Matrix float jako najwieksza bezwzgledna wartosc z jej wnetrza

        float max = fabs(m[1][1]);

        for(long i=1; i<=rows; ++i) {

            for(long j=1; j<=columns; ++j) {

                if(max < fabs(m[i][j]))
                    max = fabs(m[i][j]);
            }
        }
        return max;
    }

    Matrix transpose() {

        long n = this->rows;
        long m = this->columns;

        Matrix tmp(n,m);

        for(long i=1; i<=n; ++i) {

            for(long j=1; j<=m; ++j)

                if(i == j)
                    tmp.m[i][j] = this->m[i][j];
                else
                    tmp.m[i][j] = this->m[j][i];
        }

        return tmp;
    }

    void setName(std::string n) { name = n; }

private:

    float **m;
    long rows;
    long columns;
    std::string name;

};

class DMatrix {

public:

    ~DMatrix() { free_dmatrix(this->m, 1, rows, 1, columns); }
    DMatrix(long r = 1, long c = 1, std::string n = "") :
        m { dmatrix(1, r, 1, c) },
        rows(r),
        columns(c),
        name(n) {}

    DMatrix & operator=(const DMatrix &other){

        if(rows != other.rows || columns != other.columns) {

           std::cout << "Zle przypisanie DMatrixa!\n";
           return *this;
        }
       free_dmatrix(m, 1, rows, 1, columns);
       m = dmatrix(1, rows, 1, columns);

        for(long i = 1 ; i <= rows; i++) {
            for(long j = 1 ; j <= columns; j++)
                this->m[i][j] = other.m[i][j];
        }
       return *this;
    }

    DMatrix(const DMatrix &other){

        rows = other.rows;
        columns = other.columns;
        m = dmatrix(1, rows, 1, columns);

       for(long i = 1; i <= rows; i++) {
           for(long j = 1 ; j <= columns; j++)
               this->m[i][j] = other.m[i][j];
       }
    }

    friend std::ostream& operator<<(std::ostream& output, const DMatrix &m) {

        output << "\t\tMacierz " << m.name << "\n";

        for(long i=1; i<=m.rows; ++i) {

            for(long j=1; j<=m.columns; ++j)
                output << m(i,j) << " ";
            output << "\n";
        }
        output << "\n";

        return output;
    }

    DMatrix operator*(const DMatrix& other){ // mnozenie DMatrix przez DMatrix

        long n = this->rows;
        long m = other.columns;

        DMatrix tmp(n,m);

        for(long i = 1 ; i <= n; i++) {

            for(long j=1; j <= m; ++j) {

                tmp.m[i][j] = 0.0;
                for(long k=1; k<=this->columns; ++k)
                    tmp.m[i][j] += this->m[i][k] * other.m[k][j];
            }
        }
        return tmp;
    }

    DVector operator*(const DVector &other) { // mnozenie DMatrix przez DVector

        long n = other.getSize();

        DVector tmp(n);

        for(long i=1; i <= n; ++i)
            tmp[i] = 0.0;

        for(long i = 1 ; i <= this->rows; i++) {

            for(long j = 1; j <= this->columns; ++j)
                tmp[i] += this->m[i][j] * other[j];
        }
        return tmp;
    }

    double & operator()(long i, long j) { return m[i][j]; }
    const double & operator()(long i, long j) const { return m[i][j]; }
    double ** operator()() { return m; }

    double ptrMatrix() { // wskaznik uwarunkowania DMatrix jako najwieksza bezwzgledna wartosc z jej wnetrza

        double max = fabs(m[1][1]);

        for(long i=1; i<=rows; ++i) {

            for(long j=1; j<=columns; ++j) {

                if(max < fabs(m[i][j]))
                    max = fabs(m[i][j]);
            }
        }
        return max;
    }

    void setName(std::string n) { name = n; }

private:

    double **m;
    long rows;
    long columns;
    std::string name;

};