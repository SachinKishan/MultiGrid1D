#pragma once

#define PI 3.14159265359
#include "Eigen/Sparse"


float function2d_analytical(float x,float y)
{
    return powf(x,3) + powf(y,3);

	//return (powf(x, 2) - powf(x, 4)) * 
		//(powf(y, 4) - powf(y, 2));
}

float function2d_twicedifferentiated(float x,float y)
{
    float xcomp = 6.0f * x;
    float ycomp = 6.0f * y;

    return xcomp + ycomp;

	return 2 * ((1 - 6 * x * x) * y * y * (1 - y * y) + (1 - 6 * y * y) * x * x * (1 - x * x));
}


// Kronecker product function for sparse matrices
Eigen::SparseMatrix<float> kroneckerProduct(const Eigen::SparseMatrix<float>& A, const Eigen::SparseMatrix<float>& B) {
    int rowsA = A.rows(), colsA = A.cols();
    int rowsB = B.rows(), colsB = B.cols();

    Eigen::SparseMatrix<float> result(rowsA * rowsB, colsA * colsB);
    std::vector<Eigen::Triplet<float>> triplets;

    for (int kA = 0; kA < A.outerSize(); ++kA) {
        for (typename Eigen::SparseMatrix<float>::InnerIterator itA(A, kA); itA; ++itA) {
            for (int kB = 0; kB < B.outerSize(); ++kB) {
                for (typename Eigen::SparseMatrix<float>::InnerIterator itB(B, kB); itB; ++itB) {
                    triplets.emplace_back(itA.row() * B.rows() + itB.row(), itA.col() * B.cols() + itB.col(),
                        itA.value() * itB.value());
                }
            }
        }
    }

    result.setFromTriplets(triplets.begin(), triplets.end());
    return result;
}

void setLaplacian(Eigen::SparseMatrix<float> &A, float dx, float dy, float nx, float ny)
{
    A.setZero();
	Eigen::SparseMatrix<float> Ix = Eigen::SparseMatrix<float>(nx, nx);
	Eigen::SparseMatrix<float> Iy = Eigen::SparseMatrix<float>(ny, ny);

	Ix.setIdentity();
	Iy.setIdentity();

    // List of triplets to hold the non-zero elements
    std::vector<Eigen::Triplet<float>> tripletListx;
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
    Bhx.setZero();
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
    Bhy.setZero();
    // Set the matrix values from the triplet list
    Bhy.setFromTriplets(tripletListy.begin(), tripletListy.end());


    A= kroneckerProduct(Iy, Bhx) + kroneckerProduct(Bhy, Ix);

    //A = Bhx + Bhy;

    std::cout << A;

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

Eigen::SparseMatrix<float> createProlongationOperator(int nx)
{
    int newNx = nx * 2 - 1;
    if (newNx == 1)newNx = 2;
    Eigen::SparseMatrix<float> prolongator1D(newNx, nx);

    // Reserve enough space to avoid reallocations
    prolongator1D.reserve(Eigen::VectorXi::Constant(nx, 3)); 

    // Fill the prolongation matrix
    for (int i = 0; i < nx; i++)
    {
        // Insert the direct mapping
        prolongator1D.insert(2 * i, i) = 1.0f;  
        // Insert the interpolation if i < nx - 1
        if (i < nx - 1)
        {
            prolongator1D.insert(2 * i + 1, i) = 0.5f;   // Interpolate with the current coarse point
            prolongator1D.insert(2 * i + 1, i + 1) = 0.5f;  // Interpolate with the next coarse point
        }
    }
    return kroneckerProduct(prolongator1D,prolongator1D);

}

Eigen::VectorXf prolongate2d(Eigen::VectorXf coarseV, int nx, int ny)
{

    int newNx = nx * 2 - 1;
    if (newNx == 1)newNx = 2;
    int newNy = ny * 2 - 1;
    if (newNy == 1)newNy = 2;

    Eigen::VectorXf newMat(newNx*newNy);

    newMat.setZero();




    Eigen::SparseMatrix<float> prolongator1D(newNx, nx);

    // Reserve enough space to avoid reallocations
    //prolongator1D.reserve(Eigen::VectorXi::Constant(nx, 3)); // 2-3 non-zero elements per column

    // Fill the prolongation matrix
    for (int i = 0; i < nx; i++)
    {
        // Insert the direct mapping
        prolongator1D.insert(2 * i, i) = 1.0f;  // Direct copy of coarse grid point

        // Insert the interpolation if i < nx - 1
        if (i < nx - 1)
        {
            prolongator1D.insert(2 * i + 1, i) = 0.5f;   // Interpolate with the current coarse point
            prolongator1D.insert(2 * i + 1, i + 1) = 0.5f;  // Interpolate with the next coarse point
        }
    }

    
    std::cout << "\nPrologator\n"<<prolongator1D;

    Eigen::SparseMatrix<float> I= kroneckerProduct(prolongator1D, prolongator1D);

    Eigen::VectorXf new_v = I * coarseV;

    std::cout << "\n\ndimensions: " << I.rows() << " x " << I.cols();


    /*
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
    */
    //std::cout << "\n----new mat---------\n";

    //std::cout<<std::endl << std::endl << newMat;
    
    return new_v;

}

Eigen::VectorXf restrict2d(Eigen::VectorXf fineV, int &nx, int &ny)
{
    int newNx = (nx+1)/2  ;
    int newNy = (ny + 1) / 2 ;
    // Eigen::VectorXf newMat(newNx * newNy);
    // newMat.setZero();
    /*
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
    */
    //std::cout << "\n----restriction-----\n" << newMat;

    nx = newNx;
    ny = newNy;


    Eigen::SparseMatrix<float> I = createProlongationOperator(newNx);
    
    //std::cout << "\n\ndimensions: " << I.rows() << " x " << I.cols();
    //std::cout << "\n\n v dimensions: " << fineV.rows() << " x " << fineV.cols();

	Eigen::VectorXf new_v = I.transpose() * fineV;
    return new_v;

}

void print_as_matrix(Eigen::VectorXf v, int nx)
{
    for (int i = 0; i < nx; i++)
    {
        std::cout << "\n";

        for (int j = 0; j < nx; j++)
        {
            std::cout << v(nx*i+j) << " ";
        }
    }
}


Eigen::VectorXf multi_grid_cycle2d(Eigen::SparseMatrix<float> A, Eigen::VectorXf f, Eigen::VectorXf v, int nu1, int nu2, int nx, int ny, float dx, float dy, int level)
{
    //relaxation
    v = jacobi2D(v, f, A, nu1, 0.8f);

    //find residual
    Eigen::VectorXf residual = f - A * v;

    //restriction (to coarse)
    // the nx and ny are referenced and changed based on the restrict2d new nx and new ny
    Eigen::VectorXf coarser_residual = restrict2d(residual,nx,ny);

	Eigen::SparseMatrix<float> coarser_A(nx*nx, ny*ny);

	setLaplacian(coarser_A,dx,dy,nx,ny);


    Eigen::VectorXf coarser_v(nx*nx);
    coarser_v.setZero();
    //std::cout << "Rows of residual: " << coarser_residual.rows();
    //solve coarse grid
    if (nx > 1)
    {
        std::cout << "\nCoarsened further to an nx of ";
        std::cout << nx << "\n";
        coarser_v = multi_grid_cycle2d(coarser_A, coarser_residual, coarser_v, nu1, nu2, nx, ny, dx/2, dy/2, level+1);
    }
    else if(nx==1 || level==5)
    {
        std::cout << "max coarseness reached";
        coarser_v = jacobi2D(coarser_v, coarser_residual, coarser_A, 50, 0.8f); //solve 
    }
    
    //prolong it back
    Eigen::VectorXf prolonged_vector = prolongate2d(coarser_v, nx, ny);

    //std::cout << "\nPROLONGED: " << prolonged_vector.rows();
    //std::cout << "\ncurrent: " << v.rows();
    if (v.rows() != prolonged_vector.rows())std::cout << "\nTHE PROLONGED VECTOR AND ORIGINAL VECTOR DONT MATCH\nOriginal: " << v.rows()
        << "\nPROLONGED: " << prolonged_vector.rows() <<std::endl;
    v = v + prolonged_vector;
    //post-relaxation
    v = jacobi2D(v, f, A, nu2, 0.9f);

    return (v);

}