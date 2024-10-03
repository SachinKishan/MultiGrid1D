#pragma once

#define PI 3.14159265359
#include "Eigen/Sparse"
#include "glm.hpp"

#include "implot.h"
void set_A(Eigen::SparseMatrix<float>& A, int nx, float dx)
{
    // List of triplets to hold the non-zero elements
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(3 * nx - 2); // We have 3 diagonals: main, lower, upper

    // Fill the triplet list with values corresponding to the diagonals
    for (int i = 0; i < nx; ++i) {
        if (i > 0) {
            // Subdiagonal (below the main diagonal)
            tripletList.emplace_back(i, i - 1, -1.0 / (dx * dx));
        }

        // Main diagonal
        tripletList.emplace_back(i, i, 2.0 / (dx * dx));

        if (i < nx - 1) {
            // Superdiagonal (above the main diagonal)
            tripletList.emplace_back(i, i + 1, -1.0 / (dx * dx));
        }
    }

    // Set the matrix values from the triplet list
    A.setFromTriplets(tripletList.begin(), tripletList.end());
}


Eigen::VectorXf prolongation_vector(Eigen::VectorXf r)
{
    Eigen::VectorXf prolonged_r(r.rows() * 2 + 1);

    for (int i = 1; i < prolonged_r.rows(); i+=2)//injection step
    {
        prolonged_r(i) = r((i + 1) / 2 - 1);
    }

	prolonged_r(0) = 0.5 * prolonged_r(1);
    prolonged_r(prolonged_r.rows() - 1) = 0.5 * prolonged_r(prolonged_r.rows() - 2);

    /*
	for(int i=2;i<prolonged_r.rows()-2; i+=2)
    {
        prolonged_r(i) = 0.5 * prolonged_r(i + 1) +
            0.5 * prolonged_r(i - 1);
    }*/

    for (int i = 1; i < r.rows(); ++i)
    {
        prolonged_r(2 * i) = 0.5f * (prolonged_r(2 * i - 1) + prolonged_r(2 * i + 1));
    }



    return prolonged_r;
}

Eigen::VectorXf restriction_vector(Eigen::VectorXf r)
{

    Eigen::VectorXf restricted_vector(r.rows() / 2);
    
    for(int i=0;i<restricted_vector.rows();i++)
    {
        restricted_vector(i) = r(i*2);
    }

    return restricted_vector;
}

Eigen::SparseMatrix<float> restriction_matrix(const Eigen::SparseMatrix<float>& A) {
    int new_size = A.cols() / 2;
    Eigen::SparseMatrix<float> restricted_A(new_size, new_size);

    // List of triplets to hold the non-zero elements for restricted_A
    std::vector<Eigen::Triplet<float>> tripletList;
    tripletList.reserve(3 * new_size);  // 3 diagonals for restricted_A

    // Iterate over the rows for restricted_A
    for (int i = 0; i < new_size; ++i) {
        // Main diagonal (take every second element from the original main diagonal)
        tripletList.emplace_back(i, i, A.coeff(2 * i, 2 * i));

        // Superdiagonal (every second element from the original superdiagonal)
        if (i < new_size - 1) {
            tripletList.emplace_back(i, i + 1, A.coeff(2 * i, 2 * i + 1));
        }

        // Subdiagonal (every second element from the original subdiagonal)
        if (i > 0) {
            tripletList.emplace_back(i, i - 1, A.coeff(2 * i, 2 * i - 1));
        }
    }

    // Set the values of restricted_A from the triplet list
    restricted_A.setFromTriplets(tripletList.begin(), tripletList.end());

    return restricted_A;
}
//pass by reference, so no need for copies, no need for returns either- less memory usage

Eigen::VectorXf weighted_gauss_jacobi(Eigen::SparseMatrix<float> A, Eigen::VectorXf f, Eigen::VectorXf v, float weight, float nu = 1)
{
    Eigen::VectorXf inverse_diagonal = A.diagonal();
    Eigen::MatrixXf off_diagonal = A;
    for(int i=0;i<inverse_diagonal.rows();i++)
    {
        off_diagonal(i, i) = 0;
        inverse_diagonal(i) = 1 / inverse_diagonal(i);
    }
    
    for(int i=0;i<nu;i++)
    {
        Eigen::VectorXf residual = f - off_diagonal*v;
        v = (1 - weight) * v + weight * inverse_diagonal.cwiseProduct(residual);

    }
    
    return v;
}

Eigen::VectorXf multi_grid_cycle(Eigen::SparseMatrix<float> A, Eigen::VectorXf f, Eigen::VectorXf v, int nu1, int nu2)
{
    //relaxation
    v = weighted_gauss_jacobi(A,f,v,2.0/3.0,nu1);

    //find residual
    Eigen::VectorXf residual = f - A * v;

	//restriction (to coarse)
    Eigen::VectorXf coarser_residual = restriction_vector(residual);
    Eigen::SparseMatrix<float> coarser_A = restriction_matrix(A);
    Eigen::VectorXf coarser_v(v.rows() / 2);
    coarser_v.setZero();
    std::cout <<"Rows of residual: " << coarser_residual.rows();
    //solve coarse grid
    if(coarser_residual.rows()>1)
    {
        std::cout << "Coarsened further";
        coarser_v = multi_grid_cycle(coarser_A, coarser_residual, coarser_v, nu1, nu2);
    }
    //have a direct solve

    //prolong it back
    Eigen::VectorXf prolonged_vector = prolongation_vector(coarser_v);

    std::cout << "\nPROLONGED: " << prolonged_vector.rows();
    std::cout << "\ncurrent: " << v.rows();

	v = v + 0.25*prolonged_vector;
	//post-relaxation
    v = weighted_gauss_jacobi(A, f, v, 2.0 / 3.0,nu2);

    return (v);

}

void set_up_and_call_multigrid()
{
	
}

float function_sin(float x, float L = 1.5f, int k = 2)
{
	return sinf((float)k * PI * x / L);
}

float function_sin_diff2(float x, float L = 1.5f, int k = 2)
{
    return  pow((k * PI / L), 2) * sinf(k * PI * x / L);
}

void plot(float a, float b, std::vector<float> input, std::vector<float> output)
{
    //std::cout << "\nMax: " << max;
    ImGui::Begin("Solution Plot");
    //ImPlot::SetNextAxesLimits(a, b, -2, 2, ImGuiCond_Always);
    if (ImPlot::BeginPlot("Plot")) {
        ImPlot::PlotLine("Actual", input.data(),
            output.data(), input.size());
        ImPlot::EndPlot();
    }
    ImGui::End();
}

std::vector<float>interpolate_line(float a, float b, float step, int n)
{
    std::vector<float> output;
    //output.push_back(0);
    for (int i = 0; i < n; i++)
    {
        std::cout <<a+ i * step<<"  ";
        output.push_back(a+i * step);
    }
    return output;
}

