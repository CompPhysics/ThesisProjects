#include <mpi.h>

using namespace std;

#include "Ref/Ref.h"
#include "Domain/Domain.h"
#include "Coor/Coor.h"
#include "Vmc/Vmc.h"
#include "SlaterDet/SlaterDet.h"
#include "FuncUpDown/FuncUpDown.h"

int main(int argn, char* args[]) {

  // ************************* MPI ***********************
  int rank, size;
  MPI::Init(argn, args);
  rank = MPI::COMM_WORLD.Get_rank();
  size = MPI::COMM_WORLD.Get_size();
  system("renice +19 -u simensr");

  // Defines the form of the SlaterDet template
  typedef SlaterDet<FuncUp<CoorSpinDiff>, 
    FuncDown<CoorSpinDiff> > SlaterDeterminant;

  // ************************ Domain *********************
  Ref<Domain> domain;
  char*       initFileName = "init";
  domain      = new Domain(initFileName, rank, size);
  domain().init();

  
  // ************************* vmc ***********************
  Ref<Vmc<SlaterDeterminant> > vmc;
  vmc = new Vmc<SlaterDeterminant >(domain());
  vmc().run();

  // ************************* MPI ***********************
  MPI::Finalize();
  return 0;
}
