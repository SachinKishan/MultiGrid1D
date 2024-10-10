//#include "imgui.h"
#include <iomanip>

#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include "vendor/glad/include/glad/glad.h"
#include <GLFW/glfw3.h>
#include <iostream>
#include "Multigrid1D.h"
#include "Multigrid2D.h"
#include "spectral.h"

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow* window);

// settings
const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;



int main()
{

    // glfw: initialize and configure
    // ------------------------------
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // glfw window creation
    // --------------------
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "LearnOpenGL", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    // glad: load all OpenGL function pointers
    // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    //IMGUI INITIALISATION
    ImGui::CreateContext();

    ImGuiIO& io = ImGui::GetIO(); (void)io;
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");

    ImPlot::CreateContext();


    float a, b;//boundary values in 1D
    a = 0.0f;
    b = 1.0f;

    int nx = 8;
    float dx = (b - a) / (nx + 1);

    std::vector<float> input = interpolate_line(a + dx, b - dx, dx, nx);
    std::vector<float> output;

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            output.push_back(0);
            output[i * (nx)+j] = function2d_analytical(input[j], input[i]);
        }
    }

    //Av=f
    //laplace A
    Eigen::SparseMatrix<float> A(nx, nx);

    setLaplacian(A, dx, dx, nx, nx);
    std::cout << "\nLaplacian\n" << A;

    //f
    Eigen::VectorXf f((nx) * (nx));

    Eigen::SparseMatrix<float> F((nx) * nx, (nx)*nx);
    F.setZero();
	for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < nx; j++)
        {
            if (i == 0 || j == 0 || j == nx - 1 || i == nx - 1) //set value at boundary to known function value
            {
                f((nx)*i + j) = function2d_analytical(input[j], input[i]);
                F.coeffRef((nx)*i + j, (nx)*i + j) = 1;
            }
            else
                f((nx)*i + j) = 0; //set whatever inside as 0 ie we dont know it and the second derivative is supposed to be 0

        }
    }
	std::cout << "\nf\n";

    //forcing function
    print_as_matrix(f, nx);

    std::cout << "\nF\n";
    std::cout << F;

    Eigen::SparseMatrix<float> I = Eigen::SparseMatrix<float>(nx * nx, nx * nx);
    I.setIdentity();
    Eigen::SparseMatrix<float> L = F + (I - F) * A;

	std::cout << "\nL\n";
    std::cout << L;

    Eigen::VectorXf v(nx * nx);
    v.setZero();

    /*
    Eigen::SparseLU<Eigen::SparseMatrix<float>> solver;

    solver.compute(L);
    v = solver.solve(f);
    */

    v = multi_grid_cycle2d(L, f, v, 3, 3, nx, nx, dx, dx);


	std::cout << "\n\nFound Solution: \n";
    print_as_matrix(v, nx);

    std::cout << "\n-----Actual--------\n";
    for (int i = 0; i < nx; i++)
    {
        std::cout << "\n";

        for (int j = 0; j < nx; j++)
        {
            std::cout << output[(nx)*i + j] << " ";
        }
    }
   

	// render loop
    // -----------
    while (!glfwWindowShouldClose(window))
    {

        // render
        // ------
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);


        //after clearing buffer, declare new frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        
        ImGui::Begin("Multigrid test");
        ImGui::Text("1D Implementation");

        /*
        if(ImGui::Button("One more cycle"))
        {
            solution = multi_grid_cycle(A, f, solution, 2, 2);
        	for (int i = 0; i < solution.rows(); i++)
            {
        		solution_vector[i]=solution(i);
            }
        }
        */
        ImGui::End();
        
    	
        
        //plot(a, b, input, output);
        //plot(a, b, input, solution_vector);


        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());




        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    //end imgui functions
    ImPlot::DestroyContext();
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();


    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    glfwTerminate();
    
    return 0;
}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and 
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}