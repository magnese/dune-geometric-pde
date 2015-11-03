#define GRIDDIM ALBERTA_DIM-1
#define WORLDDIM ALBERTA_DIM

#define ALBERTAGRID // the mesh file must be oriented consistently since it does not perform any check
//#define ALUGRID_SIMPLEX

#define HAVE_UMFPACK 1

#define POLORDER 1

#include "config.h"
#include <dune/common/timer.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/gmshwriter.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/geometrygrid/grid.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/io/parameter.hh>

#include <string>
#include <iostream>
#include <vector>

#include "vertexfunction.hh"
#include "femschemeinterface.hh"
#include "computeinterface.hh"

int main(int argc,char** argv)
{
  try
  {
    // start timer
    Dune::Timer timer(false);
    timer.start();

    // init
    Dune::Fem::MPIManager::initialize(argc,argv);
    Dune::Fem::Parameter::append(argc,argv);
    Dune::Fem::Parameter::append((argc<2)?("/home/ma/m/ma2413/dune-repo/dune-geometric-pde/src/parameter"):(argv[1]));

    // load host grid
    typedef Dune::GridSelector::GridType HostGridType;
    const std::string fileName(Dune::Fem::Parameter::getValue<std::string>("FileName","mesh.msh"));
    Dune::GridFactory<HostGridType> hostGridFactory;
    std::vector<int> boundaryIDs(0);
    std::vector<int> elementsIDs(0);
    Dune::GmshReader<HostGridType>::read(hostGridFactory,fileName,boundaryIDs,elementsIDs);
    HostGridType* hostGrid(hostGridFactory.createGrid());

    // create grid
    typedef Dune::GeometryGrid<HostGridType,Dune::Fem::VertexFunction<HostGridType>> GridType;
    GridType grid(hostGrid);

    // load parameter
    const bool useMeanCurvatureFlow(Dune::Fem::Parameter::getValue<bool>("UseMeanCurvatureFlow",0));
    const std::string fileNameFinalMesh(Dune::Fem::Parameter::getValue<std::string>("FileNameFinalMesh",""));
    if(useMeanCurvatureFlow)
      std::cout<<"Problem type: mean curvature flow."<<std::endl;
    else
      std::cout<<"Problem type: surface diffusion."<<std::endl;

    // compute solution
    typedef Dune::Fem::FemSchemeInterface<GridType> FemSchemeType;
    FemSchemeType femScheme(grid,useMeanCurvatureFlow);
    Dune::Fem::computeInterface<FemSchemeType>(femScheme);

    // dump final mesh as msh
    if(fileNameFinalMesh!="")
    {
      Dune::GmshWriter<typename GridType::LeafGridView> gmshWriter(grid.leafGridView());
      //gmshWriter.setPrecision(15);
      gmshWriter.write(fileNameFinalMesh,elementsIDs);
      std::cout<<std::endl<<"Final mesh dumped into "<<fileNameFinalMesh<<"."<<std::endl;
    }

    // output total running time
    timer.stop();
    std::cout<<std::endl<<"Total running time: "<<timer.elapsed()<<" seconds."<<std::endl;

    return 0;
  }

  catch(Dune::Exception &e)
  {
    std::cerr<<"Dune reported error: "<<e<<std::endl;
  }

  catch(...)
  {
    std::cerr<<"Unknown exception thrown!"<<std::endl;
  }

}
