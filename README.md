Newton's Method implementation for solving convex quadratic programs via the barrier method.

We want to minimize a function that can be written as 1/2*x^T*Q*x+c^T*x subject to the conditions A*x=b, x≥0, where Q and A are given matrices and c is a given vector.

Requirements: Matlab (for base implementation) and either Yalmip and SDP3 (for testing). It can easily be modified to use a test from the optimization toolbox.

The file OptimizationTest.m generates a random constraint problem and runs the implementation against a standard solver. The test takes a number of user inputs: n gives the number of dimensions, m the number of equality constraints, and condnum gives the condition number on the randomly generated matrix. S_0 and mu control how much gain each iteration must produce, thus controling the number of iterations required. Note that there is a tradeoff with the time required to run each iteration. Finally, epsilon and maxiter give stopping conditions.
OptimizationTest saves a file OptimizationOutput.txt detailing the optimal solution produced by the algorithm, the optimal function value, and the number of iterations and newton steps required to reach the solution. It also gives the optimal solution given by Yalmip and SDP3 as a comparison.

The implementation itself is run through the function Barrier.m.

Barrier takes the function to be optimized F, its gradient GradF, its hessian HessF, along with a feasible initial point x_0, the constraint matrix A, several constants controling speed/accuracy tradeoffs and final approximation accuracy (s_0, mu, epsilon, alpha, beta, and maxiter), Q and c from the problem setup, and the filepath to which output should be saved (fid).

s_0 and mu control the accuracy required of sub steps in this implementation. High s_0 and mu will decrease the number of iterations required, but the iterations will be slower. The most efficient s_0 and mu depend on the specific problem being addressed. Epsilon is the accuracy requirement for the problem as a whole. The program returns a function value that is within epsilon of the optimal value. Epsilon is required to be positive. Relatively large epsilon makes the program run quickly but gives a poor approximation, while small epsilon is accurate but slow. Alpha and Beta also determine a speed/accuracy trade-off, here in the calculation of newton steps. Alpha is between 0 and .5 and beta is between 0 and 1. Higher alpha and beta make the program require fewer newton steps, but make the steps take longer to calculate. Finally, maxiter is the maximum number of iterations the program will take before giving up. A high maxiter may solve problems that a low maxiter rejects, but it takes far longer to detect an error.

The program returns Inneriter (the number of newton steps taken), Solutioniter (the number of constrained newton iterations taken), Solutionval (the optimal function value to within epsilon), and Solutionx (the vector which produces said function value).
