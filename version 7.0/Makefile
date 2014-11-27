EXED = build/
SRCD = src/
PETSC_DIR = /home/ninoy/libs/petsc-3.4.3
PETSC_ARCH = arch-linux2-cxx-debug
PAPI_DIR = /usr/local

FC = mpif90.mpich2

FFLAGS = -g -cpp -fPIC -fopenmp -mcmodel=large

INC_PATH = -I$(PETSC_DIR)/include/ -I$(PETSC_DIR)/$(PETSC_ARCH)/include -I$(PAPI_DIR)/include/

LIB_PATH = -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -L$(PAPI_DIR)/lib

LIBS = -lpapi -lpetsc -Wl,-rpath=$(PAPI_DIR)/lib -Wl,-rpath=$(PETSC_DIR)/$(PETSC_ARCH)/lib
#-lHYPRE -lspai 

EXECTBL=GravitySolver

SRCS=MD_Definition.f90 MD_Parameter.f90 MD_Quantity.f90 MD_GeometricQuantity.f90 MD_PETScQuantity.f90 \
MD_IO.f90 MD_Helper.f90 DomainDecomposition.f90 MD_CalculateCordinate.f90 MD_BoundaryCondition.f90 \
IC_Vaccum.f90 IC_UniformDensity.f90 IC_SinusoidalDensity.f90 IC_Poly6Density.f90 \
IC_CentredSphere.f90 IC_OffCentredGaussian.f90 IC_GaussianDensity.f90 \
IC_CondensedSphere.f90 IC_OffCentredCondensed.f90 IC_Ploytrop1.f90 IC_Ploytrop5.f90 \
ComputeCartesianMapping.f90 ComputeTransverseGradient.f90 ComputeAvgDensity.f90 \
ComputeStencil2.f90 ComputeStencil4.f90 ComputeMatrix.f90 ComputeRHS.f90 \
Test_Indexing.f90 Test_BoundaryValue.f90 Test_Accuracy.f90 Test_Performance.f90 \
Test_SolverConvergence.f90 RunTest.f90 \
Init.f90 Solver.f90 Output_1D.f90 Output_3D.f90 Finalize.f90 \
ParallelGravitySolver.f90

OBJS=$(SRCS:.f90=.o)

BINOBJECTS=$(patsubst %.o,$(EXED)%.o,$(OBJS))

all : makedir $(EXECTBL) move

$(EXECTBL): $(OBJS)
	@echo ''
	@echo '============ ...building the executable ============'
	@echo ''
	$(FC) $(FFLAGS) $(INC_PATH) $(BINOBJECTS) -o $(EXED)$@ $(LIB_PATH) $(LIBS)
	@echo ''
	@echo '============ ...Done! ============'

$(OBJS):%.o : $(SRCD)%.f90
	@echo ''
	@echo '============ ...building' $*.o '============'
	@echo ''
	$(FC) $(FFLAGS) $(INC_PATH) -c $< -o $(EXED)$@

makedir:
	mkdir -p build data

move:
	mv *.mod $(EXED)

clean:
	rm -f $(BINOBJECTS) $(EXED)$(EXECTBL) $(EXED)*.mod

flush: 
	rm -f data/*.dat data/*.vtk

flushlog:
	rm -f data/*.log
