EXED = build/
SRCD = src/
PETScD = /home/ninoy/libs/petsc-3.4.3

FC = mpif90.mpich2

FFLAGS = -g -cpp -fPIC -fopenmp -mcmodel=large 
#FFLAGS = -g -cpp -fPIC -mcmodel=large -pg

INC_PATH = -I$(PETScD)/include/ -I$(PETScD)/arch-linux2-cxx-debug/include

LIB_PATH = -L$(PETScD)/arch-linux2-cxx-debug/lib

LIBS = -lpetsc

EXECTBL=GravitySolver


SRCS=MD_Definition.f90 MD_Parameter.f90 MD_Quantity.f90 MD_GeometricQuantity.f90 \
MD_CalculateCordinate.f90 MD_CordinateTransform.f90 \
MD_IO.f90 MD_Helper.f90 MD_BoundaryCondition.f90 \
Test_Indexing.f90 Test_BoundaryValue.f90 Test_Accuracy.f90 Test_SolverConvergence.f90 RunTest.f90 \
IC_Vaccum.f90 IC_UniformDensity.f90 IC_SinusoidalDensity.f90 IC_Poly6Density.f90 \
IC_CentredSphere.f90 IC_OffCentredSphere.f90 IC_GaussianDensity.f90 Init.f90 \
ComputeTransverseGradient.f90 CalculateStencil.f90 ComputeMatrix.f90 \
CalculateAvgDensity.f90 ComputeRHS.f90 Solver.f90 ComputePointValue.f90 \
Output_1D.f90 Output_3D.f90 OutputError.f90 Finalize.f90 GravitySolver.f90 

OBJS=$(SRCS:.f90=.o)

BINOBJECTS=$(patsubst %.o,$(EXED)%.o,$(OBJS))

all : $(EXECTBL) move

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
	mkdir -p build
	mkdir -p data
	$(FC) $(FFLAGS) $(INC_PATH) -c $< -o $(EXED)$@

move:
	mv *.mod $(EXED)

clean:
	rm -f $(BINOBJECTS) $(EXED)$(EXECTBL) $(EXED)*.mod

flush: 
	rm -f data/*.dat data/*.vtk
