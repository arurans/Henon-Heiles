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

int main()
{
    vec y0 = create_init_cond(H_0);
    //compute_hamiltonians(t_0, t_end, y0, h);
    //compute_poincare_maps(t_0, t_end, y0, h);
    compute_both(t_0, t_end, y0, h);

	return 0;
}