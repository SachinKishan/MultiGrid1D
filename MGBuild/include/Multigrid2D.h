#pragma once

#define PI 3.14159265359
#include "Eigen/Sparse"
#include "glm.hpp"

#include "implot.h"

float function2d_analytical(float x,float y)
{
	return (powf(x, 2) - powf(x, 4)) * 
		(powf(y, 4) - powf(y, 2));
}

float function2d_twicedifferentiated(float x,float y)
{
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
    Bhx.setFromTriplets(tripletListy.begin(), tripletListy.end());

    A = Bhx + Bhy;

    //std::cout << A;

}

void jacobi2D()
{
}

/**
 * \brief 
 * \param coarseV 2n x 1 vector or n x 2
 * \param nx 
 * \param ny 
 */
void prolongate2d(Eigen::MatrixXf coarseV, int nx, int ny)
{

    int newNx = nx * 2 - 1;
    int newNy = ny * 2 - 1;
    Eigen::MatrixXf newMat(newNx*newNy,2);

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
            newMat.row(newNx * icoord + jcoord) =coarseV.row(nx*i+j);
            if (i < nx - 1)
                newMat.row(newNx * (icoord + 1) + jcoord) = 0.5 * (coarseV.row(nx*i+j) + coarseV.row(nx * (i+1) + j));
            if (j < ny - 1)
                newMat.row(newNx * icoord + jcoord + 1) = 0.5 * (coarseV.row(nx * i + j) + coarseV.row(nx * i + j+1));
            if (i < nx - 1 && j < ny - 1)
                newMat.row(newNx * (icoord + 1) + jcoord + 1) = 0.25 * (coarseV.row(nx * i + j) + coarseV.row(nx * (i + 1) + j) + coarseV.row(nx * i + j + 1)
                     + coarseV.row(nx * (i + 1) + j + 1));
	    }
    }
    std::cout << "\n----new mat---------\n";

    std::cout<<std::endl << std::endl << newMat;


}

void restrict2d()
{
	
}