#include "compute.h"

/**
 * @brief 
 * Compute the hamiltonians for all the implemented methods,
 * and save them in a file
 * 
 * @param t_0 start time
 * @param t_end end time
 * @param y0 initial condition for system
 * @param h length of timestep
 */
void compute_hamiltonians(const double& t_0, const double& t_end, const vec& y0, const double& h)
{
    //Time array
    vec T = create_T(t_0, t_end, h);

    //Matrix containing the time in the first column and hamiltonians of all methods in the second
    //Should have n rows (number of time steps) and the number of methods + 1 columns
    mat H = zeros(T.n_elem, 5);
    H.unsafe_col(0) = T;

    #pragma omp parallel
    {
        #pragma omp single nowait
        {
            // Compute the hamiltonian of Kutta's method
            #pragma omp task
            H.unsafe_col(1) = hamiltonian(kuttas_method(t_0, t_end, y0, h));

            // Compute the hamiltonian of Shampine-Bogacki
            #pragma omp task
            H.unsafe_col(2) = hamiltonian(shampine_bogacki(t_0, t_end, y0, h));

            // Compute the hamiltonian of Kahans method
            #pragma omp task
            H.unsafe_col(3) = hamiltonian(kahans(t_0, t_end, y0, h));

            // Compute the hamiltonian of Störmer-Verlet
            #pragma omp task
            H.unsafe_col(4)  = hamiltonian(stormer_verlet(t_0, t_end, y0, h));
        }
    }
    #pragma omp taskwait

    //Uncomment the line below if you actaully want the output
    // H.save(hamiltonians_file + decimal_to_string(h), csv_ascii);
}

/**
 * @brief 
 * Compute the Poincaré maps of all the implemented methods,
 * and save each of them in a seperate file, because of
 * potential difference in array sizes
 * 
 * @param t_0 start time
 * @param t_end end time
 * @param y0 intial condition
 * @param h step size
 */
void compute_poincare_maps(const double& t_0, const double& t_end, const vec& y0, const double& h)
{
    //Declare the matrices to store results in
    mat P_rk, P_sb, P_kahans, P_sv;

    #pragma omp parallel
    {
        #pragma omp single nowait
        {
            // Find the Poincaré map of Kutta's method
            #pragma omp task
            P_rk = poincare(kuttas_method(t_0, t_end, y0, h));

            // Find the Poincaré map of Shampine-Bogacki
            #pragma omp task
            P_sb = poincare(shampine_bogacki(t_0, t_end, y0, h));

            // Find the Poincaré map of Kahans method
            #pragma omp task
            P_kahans = poincare(kahans(t_0, t_end, y0, h));

            // Find the Poincaré map of Störmer-Verlet
            #pragma omp task
            P_sv = poincare(stormer_verlet(t_0, t_end, y0, h));
        }
    }
    #pragma omp taskwait

    //Uncomment the line below if you actaully want the output
    // P_rk.save(poincare_file_rk4 + decimal_to_string(h), csv_ascii);

    //Uncomment the line below if you actaully want the output
    // P_sb.save(poincare_file_sb + decimal_to_string(h), csv_ascii);

    //Uncomment the line below if you actaully want the output
    // P_kahans.save(poincare_file_kahans + decimal_to_string(h), csv_ascii);

    //Uncomment the line below if you actaully want the output
    // P_sv.save(poincare_file_sv + decimal_to_string(h), csv_ascii);
}

/**
 * @brief 
 * Compute both the hamiltonians and Poincaré maps of
 * all the implemented methods
 * 
 * @param t_0 start time
 * @param t_end end time
 * @param y0 initial condition for system
 * @param h length of timestep
 */
void compute_both(const double& t_0, const double& t_end, const vec& y0, const double& h)
{
    //Matrices to store 
    mat Y_rk, Y_sb, Y_kahans, Y_sv;
    mat P_rk, P_sb, P_kahans, P_sv;

    //Compute the Hénon-Heiles system for each of the implemented method
    //in parallel first, as each of them will be used twice in the 
    //subsequent computation part
    #pragma omp parallel
    {
        #pragma omp single nowait
        {
            // Compute the Hénon-Heiles system with Kutta's method
            #pragma omp task
            Y_rk = kuttas_method(t_0, t_end, y0, h);

            // Compute the Hénon-Heiles system with Shampine-Bogacki
            #pragma omp task
            Y_sb = shampine_bogacki(t_0, t_end, y0, h);

            // Compute the Hénon-Heiles system with Kahans method
            #pragma omp task
            Y_kahans = kahans(t_0, t_end, y0, h);

            // Compute the Hénon-Heiles system with Störmer-Verlet
            #pragma omp task
            Y_sv = stormer_verlet(t_0, t_end, y0, h);
        }
    }
    #pragma omp taskwait

    //Time array
    vec T = create_T(t_0, t_end, h);

    //Matrix containing the time in the first column and hamiltonians of all methods in the second
    //Should have n rows (number of time steps) and the number of methods + 1 columns
    mat H = zeros(T.n_elem, 5);

    //Add the time vector to our matrix
    H.unsafe_col(0) = T;

    //Compute the hamiltonians and the Poincaré maps for each implemented method in parallel
    #pragma omp parallel
    {
        #pragma omp single nowait
        {
            // Compute the hamiltonian of Kutta's method
            #pragma omp task
            H.unsafe_col(1) = hamiltonian(Y_rk);

            // Compute the hamiltonian of Shampine-Bogacki
            #pragma omp task
            H.unsafe_col(2) = hamiltonian(Y_sb);

            // Compute the hamiltonian of Kahans method
            #pragma omp task
            H.unsafe_col(3) = hamiltonian(Y_kahans);

            // Compute the hamiltonian of Störmer-Verlet
            #pragma omp task
            H.unsafe_col(4) = hamiltonian(Y_sv);

            // Find the Poincaré map of Kutta's method
            #pragma omp task
            P_rk = poincare(Y_rk);

            // Find the Poincaré map of Shampine-Bogacki
            #pragma omp task
            P_sb = poincare(Y_sb);

            // Find the Poincaré map of Kahans method
            #pragma omp task
            P_kahans = poincare(Y_kahans);

            // Find the Poincaré map of Störmer-Verlet
            #pragma omp task
            P_sv = poincare(Y_sv);
        }
    }
    #pragma omp taskwait

    //Uncomment the line(s) below if you actaully want the output
    // H.save(hamiltonians_file + decimal_to_string(h), csv_ascii);

    // P_rk.save(poincare_file_rk4 + decimal_to_string(h), csv_ascii);

    // P_sb.save(poincare_file_sb + decimal_to_string(h), csv_ascii);

    // P_kahans.save(poincare_file_kahans + decimal_to_string(h), csv_ascii);

    // P_sv.save(poincare_file_sv + decimal_to_string(h), csv_ascii);
}