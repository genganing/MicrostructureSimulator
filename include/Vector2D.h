#ifndef VECTOR2D_H
#define VECTOR2D_H

#include <iostream>
#include "vector"

template<typename _Ty>
class Vector2D
{
    public:

    int Nx;
    int Ny;

    std::vector<_Ty> vector2D;

    Vector2D(){};
    Vector2D(int const &ni, int const &nj)
    {
        this->Nx = ni;
        this->Ny = nj;
        vector2D.resize(Nx*Ny);
    }
    
    void Resize(int const &ni, int const &nj)
    {
        this->Nx = ni;
        this->Ny = nj;
        vector2D.resize(Nx * Ny);
    }

    _Ty const & operator()(int const &i, int const &j) const
    {
        return (vector2D[j*this->Nx + i]);
    }

    _Ty & operator()(int const &i, int const &j) 
    {
        return (vector2D[j*this->Nx + i]);
    }

    _Ty const * operator[](int const &j) const
    {
        return &(vector2D[j*this->Nx]);
    }

    _Ty * operator[](int const &j) 
    {
        return &(vector2D[j*this->Nx]);
    }
    ~Vector2D() {
    };
};
#endif