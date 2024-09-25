#pragma once

#define PI 3.14159265359
#include "Eigen/Sparse"
#include "glm.hpp"

#include "implot.h"

float function2d_analytical(float x,float y)
{
    return powf(x,3) + powf(y,2);

	//return (powf(x, 2) - powf(x, 4)) * 
		//(powf(y, 4) - powf(y, 2));
}

float function2d_twicedifferentiated(float x,float y)
{
    float xcomp = 6.0f * x;
    float ycomp = 2.0f;

    return xcomp + ycomp;

	return 2 * ((1 - 6 * x * x) * y * y * (1 - y * y) + (1 - 6 * y * y) * x * x * (1 - x * x));
}

void setLaplacian(Eigen::SparseMatrix<float> &A, float dx, float dy, float nx, float ny)
{

	Eigen::SparseMatrix<float> Ix = Eigen::SparseMatrix<float>(nx, nx);
	Eigen::SparseMatrix<float> Iy = Eigen::SparseMatrix<float>(ny, ny);

	Ix.setIdentity();
	Iy.setIdentity();

    // List of triplets to hold the non-zero elements
    std::vector<Eigen::Triplet<double>> tripletListx;
    tripletListx.reserve(3 * nx - 2); // We have 3 diagonals: main, lower, upper

    // Fill the triplet list with values corresponding to the diagonals
    for (int i = 0; i < nx; ++i) {
        if (i > 0) {
            // Subdiagonal (below the main diagonal)
            tripletListx.emplace_back(i, i - 1, -1.0 / (dx * dx));
        }

        // Main diagonal
        tripletListx.emplace_back(i, i, 2.0 / (dx * dx));

        if (i < nx - 1) {
            // Superdiagonal (above the main diagonal)
            tripletListx.emplace_back(i, i + 1, -1.0 / (dx * dx));
        }
    }

    Eigen::SparseMatrix<float> Bhx(nx,ny);
    // Set the matrix values from the triplet list
    Bhx.setFromTriplets(tripletListx.begin(), tripletListx.end());

    // List of triplets to hold the non-zero elements
    std::vector<Eigen::Triplet<float>> tripletListy;
    tripletListy.reserve(3 * nx - 2); // We have 3 diagonals: main, lower, upper

    // Fill the triplet list with values corresponding to the diagonals
    for (int i = 0; i < nx; ++i) {
        if (i > 0) {
            // Subdiagonal (below the main diagonal)
            tripletListy.emplace_back(i, i - 1, -1.0 / (dx * dx));
        }

        // Main diagonal
        tripletListy.emplace_back(i, i, 2.0 / (dx * dx));

        if (i < nx - 1) {
            // Superdiagonal (above the main diagonal)
            tripletListy.emplace_back(i, i + 1, -1.0 / (dx * dx));
        }
    }

    Eigen::SparseMatrix<float> Bhy(nx, ny);
    // Set the matrix values from the triplet list
    Bhy.setFromTriplets(tripletListy.begin(), tripletListy.end());

    A = Bhx + Bhy;

    //std::cout << A;

}

Eigen::VectorXf jacobi2D(Eigen::VectorXf v, Eigen::VectorXf f,Eigen::SparseMatrix<float>& A,int n, float w)
{
    Eigen::VectorXf inverse_diagonal = A.diagonal();
    Eigen::MatrixXf off_diagonal = A;
    for (int i = 0; i < inverse_diagonal.rows(); i++)
    {
        off_diagonal(i, i) = 0;
        inverse_diagonal(i) = 1 / inverse_diagonal(i);
    }

    for (int i = 0; i < n; i++)
    {
    	Eigen::VectorXf residual = f - off_diagonal * v;
    	v = (1 - w) * v + w * inverse_diagonal.cwiseProduct(residual);
    }
    return v;
}

Eigen::VectorXf prolongate2d(Eigen::VectorXf coarseV, int nx, int ny)
{

    int newNx = nx * 2 - 1;
    if (newNx == 1)newNx = 2;
    int newNy = ny * 2 - 1;
    if (newNy == 1)newNy = 2;

    Eigen::VectorXf newMat(newNx*newNy);

    newMat.setZero();

    std::cout << "\n----Original mat---------\n";
    std::cout << coarseV;
    std::cout << "\n----\n";
    for(int i=0;i<nx;i++)
    {
	    for(int j=0;j<ny;j++)
	    {

            int icoord = i * 2;
            int jcoord = j * 2;

            std::cout << std::endl << i << " " << j<<" "<< newNx * icoord + jcoord;
            newMat(newNx * icoord + jcoord) =coarseV(nx*i+j);
            if (i < nx - 1)
                newMat(newNx * (icoord + 1) + jcoord) = 0.5f * (coarseV(nx*i+j) + coarseV(nx * (i+1) + j));
            if (j < ny - 1)
                newMat(newNx * icoord + jcoord + 1) = 0.5f * (coarseV(nx * i + j) + coarseV(nx * i + j+1));
            if (i < nx - 1 && j < ny - 1)
                newMat(newNx * (icoord + 1) + jcoord + 1) = 0.25f * (coarseV(nx * i + j) + coarseV(nx * (i + 1) + j) + coarseV(nx * i + j + 1)
                     + coarseV(nx * (i + 1) + j + 1));
	    }
    }
    std::cout << "\n----new mat---------\n";

    std::cout<<std::endl << std::endl << newMat;
    return newMat;

}

Eigen::VectorXf restrict2d(Eigen::VectorXf fineV, int &nx, int &ny)
{
    int newNx = (nx+1)/2  ;
    int newNy = (ny + 1) / 2 ;
    Eigen::VectorXf newMat(newNx * newNy);
    newMat.setZero();

    for(int i=0;i<newNx;i++)
    {
	    for(int j=0;j<newNy;j++)
	    {
            int icoord = i * 2;
            int jcoord = j * 2;
            float a=0, b=0, c=0, d=0;//4 contributing values
            
            if(i>0 && j>0)
            	a = fineV((nx) * (icoord - 1) + jcoord - 1);
            if (i > 0)
                b = fineV((nx) * (icoord - 1) + jcoord);
            if(j>0)
                c= fineV((nx) * icoord + jcoord - 1);
            d = fineV((nx) * icoord + jcoord);

            newMat(newNx * i + j) = 1.0f / 16.0f * (a + 4 * b + 4 * c + 16 * d);
	    }
    }
    std::cout << "\n----restriction-----\n" << newMat;

    nx = newNx;
    ny = newNy;

    return newMat;

}

void print_as_matrix()
{
	
}


Eigen::VectorXf multi_grid_cycle2d(Eigen::SparseMatrix<float> A, Eigen::VectorXf f, Eigen::VectorXf v, int nu1, int nu2, int nx, int ny, float dx, float dy)
{
    //relaxation
    v = jacobi2D(v, f, A, nu1, 0.8f);

    //find residual
    Eigen::VectorXf residual = f - A * v;

    //restriction (to coarse)

    Eigen::VectorXf coarser_residual = restrict2d(residual,nx,ny);

	Eigen::SparseMatrix<float> coarser_A(nx*nx, ny*ny);

	setLaplacian(coarser_A,dx,dy,nx*nx,ny*ny);


    Eigen::VectorXf coarser_v(nx*nx);
    coarser_v.setZero();
    //std::cout << "Rows of residual: " << coarser_residual.rows();
    //solve coarse grid
    if (nx > 1)
    {
        std::cout << "Coarsened further";
        std::cout << "\n" << nx << "\n";
        coarser_v = multi_grid_cycle2d(coarser_A, coarser_residual, coarser_v, nu1, nu2, nx, ny, dx/2, dy/2);
    }
    
    //prolong it back
    Eigen::VectorXf prolonged_vector = prolongate2d(coarser_v, nx, ny);

    //std::cout << "\nPROLONGED: " << prolonged_vector.rows();
    //std::cout << "\ncurrent: " << v.rows();

    v = v + prolonged_vector;
    //post-relaxation
    v = jacobi2D(v, f, A, nu2, 0.8f);

    return (v);

}