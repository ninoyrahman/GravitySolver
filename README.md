4th order poisson solver in mapped coordinate (colella 2011). This code solves poisson equation for gravity. 
First the poisson equation is solved in the cartesian computational domain then mapped to physically appropiate domain. 
The code contains both 2nd and 4th order colella 2011 stencil for elliptic equation. It has been tested for several test problems. 
Mapping is not tested yet. Parallelization is done by MPI and OpenMP. 
