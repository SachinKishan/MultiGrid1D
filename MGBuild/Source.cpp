//#include "imgui.h"
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

    int max_levels = 3; //the finest level it will go, also the level it starts from
    int nx = pow(2.0, max_levels) - 1;

    //number of cycles
    int mu = 2;

    //number of pre and post relaxations
    int pre_nu = 1;
    int post_nu = 1;

    nx = 5;
    float dx = (b - a) / (nx + 1.0f);

    //std::cout << nx << std::endl;

    std::vector<float> input= interpolate_line(a + dx, b - dx, dx, nx);
    std::vector<float> output;

    for(int i=0;i<nx;i++)
    {
        for(int j=0;j<nx;j++)
        output.push_back(function2d_analytical(input[j], input[i]));
    }

    
    
    

    
    //Av=f
    //laplace A
    Eigen::SparseMatrix<float> A(nx, nx);

    setLaplacian(A, dx, dx, nx, nx);

    //v
    //nx = 2;
    Eigen::VectorXf v(nx * nx);
    v.setZero();
    //v(0) = 1;
    //v(4) = 1;
    //v(8) = 1;

    //print_as_matrix(v, nx);

    //v = prolongate2d(v, nx, nx);
    //std::cout << "\n-------\n";
  //  print_as_matrix(v,nx*2-1);


    //v=restrict2d(v, nx, nx);
    //std::cout << "\n-------\n";
    //print_as_matrix(v,nx);
    
    
    //f
    Eigen::VectorXf f(nx*nx);
    for(int i=0;i<nx;i++)
    {
        for (int j = 0; j < nx; j++)
        {
            f(nx*i+j) = function2d_twicedifferentiated(input[j], input[i]);
            //f(nx*i+j) = function2d_analytical(input[j], input[i]);
        }
    }
    //std::cout << std::endl<<f;


    /*
    Eigen::VectorXf solution = multi_grid_cycle(A, f, v, 2, 2);
    //solution = v;
    std::cout << solution;

    std::vector<float> solution_vector;
    for(int i=0;i<solution.rows();i++)
    {
        solution_vector.push_back(solution(i));
    }
    */

    //Eigen::VectorXf solution = 
    
	v=multi_grid_cycle2d(A,f,v,2,2,nx,nx,dx,dx,1);


    std::cout << "\n-----Real--------\n";
    for (int i = 0; i < nx; i++)
    {
        std::cout << "\n";

        for (int j = 0; j < nx; j++)
        {
            std::cout << output[nx*i+j] << " ";
        }
    }
    std::cout << "\n-----Found--------\n";

    print_as_matrix(v,nx);
    
	// render loop
    // -----------
    while (!glfwWindowShouldClose(window))
    {
        // input
        // -----
        processInput(window);

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

    glfwTerminate();
    
    return 0;
}

// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------
void processInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and 
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}