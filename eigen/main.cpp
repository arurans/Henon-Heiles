#include "./src/problems/compute.h"
#include "./src/constants.h"


/**
 * I really am not interested in the output,
 * and am only aiming for fast computation (^:
 * 
 * If you for some (crazy) reason are interested
 * in the output, uncomment the proper lines of 
 * code in ./src/problems/compute.cpp
 */

/**
 * To not unnecessarily store every value used for computation,
 * there is a variable in ./src/storage_info.cpp called 
 * SKIP_STORAGE that lets the functions store only every
 * "SKIP_STORAGE-th" value.
 * 
 * This allows for computation for longer time intervals too,
 * without needing too much RAM, in addition to speeding up
 * the process of storing the data to a CSV file.
 */

#include <iostream> 

int main()
{
    Array<double, 4, 1> y0 = create_init_cond(H_0);
    //compute_hamiltonians(t_0, t_end, y0, h);
    //compute_poincare_maps(t_0, t_end, y0, h);
    compute_both(t_0, t_end, y0, h);

    return 0;
}