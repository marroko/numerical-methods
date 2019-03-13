#pragma once

#include "/opt/NR/numerical_recipes.c/nrutil.h"
#include "/opt/NR/numerical_recipes.c/nrutil.c"
#include "/opt/NR/numerical_recipes.c/ludcmp.c"
#include "/opt/NR/numerical_recipes.c/lubksb.c"

#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>

class Matrix {

public:

    ~Matrix() { free_matrix(this->m, 1, rows, 1, columns); }
    Matrix(int r = 1, int c = 1, std::string n = "") :
                                                m { matrix(1, r, 1, c) },
                                                rows(r),
                                                columns(c),
                                                name(n) {}

    Matrix & operator=(const Matrix &other){

       rows = other.rows;
       columns = other.columns;
       m = matrix(1, rows, 1, columns);

           for(int i = 0 ; i < rows; i++) {
               for(int j = 0 ; j < columns; j++)
                   this->m[i][j] = other.m[i][j];
           }
       return *this;
    }

    Matrix & operator=(Matrix &other){

        if(this == &other)
            return *this;

       rows = other.rows;
       columns = other.columns;

           for(int i = 0 ; i < rows; i++) {
               for(int j = 0 ; j < columns; j++)
                   this->m[i][j] = other.m[i][j];
           }
       return *this;
    }


    friend std::ostream& operator<<(std::ostream& output, Matrix &m) {

        output << "\t\tMacierz " << m.name << "\n";

        for(int i=1; i<=m.rows; ++i) {

            for(int j=1; j<=m.columns; ++j)
                output << m(i,j) << " ";
            output << "\n";
        }
        output << "\n";

        return output;
    }

    Matrix * operator*(Matrix& other){

        int n = this->rows;
        int m = other.columns;

        Matrix *tmp = new Matrix(n,m);

        for(int i = 1 ; i <= n; i++) {

            for(int j=1; j <= m; ++j) {

                tmp->m[i][j] = 0.0;
                for(int k=1; k<=this->columns; ++k)
                    tmp->m[i][j] += this->m[i][k] * other.m[k][j];
            }
        }
        return tmp;
    }

    float & operator()(const int i, const int j) { return m[i][j]; }
    float ** operator()() { return m; }

    float ptrMatrix() {

        float max = fabs(m[1][1]);

        for(int i=1; i<=rows; ++i) {

            for(int j=1; j<=columns; ++j) {

                if(max < fabs(m[i][j]))
                    max = fabs(m[i][j]);
            }
        }
        return max;
    }

    void setName(std::string n) { name = n; }

private:

    float **m;
    int rows;
    int columns;
    std::string name;

};