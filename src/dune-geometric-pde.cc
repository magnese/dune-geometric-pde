#define POLORDER 1

#include "config.h"
#include <dune/common/timer.hh>
#include <dune/grid/albertagrid.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/gmshwriter.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/geometrygrid/grid.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/io/io.hh>
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
    Dune::Fem::Parameter::append(argc<2?(static_cast<std::string>(SOURCEDIR)+"/src/parameter"):argv[1]);

    // load host grid
    typedef Dune::AlbertaGrid<GRIDDIM,WORLDDIM> HostGridType;
    const std::string fileName(static_cast<std::string>(MSHFILESDIR)+Dune::Fem::Parameter::getValue<std::string>("FileName","mesh.msh"));
    Dune::GridFactory<HostGridType> hostGridFactory;
    std::vector<int> boundaryIDs(0);
    std::vector<int> elementsIDs(0);
    Dune::GmshReader<HostGridType>::read(hostGridFactory,fileName,boundaryIDs,elementsIDs);
    HostGridType* hostGrid(hostGridFactory.createGrid());

    // create grid
    typedef Dune::GeometryGrid<HostGridType,Dune::Fem::VertexFunction<HostGridType>> GridType;
    GridType grid(hostGrid);

    // load problem type
    const bool useMeanCurvatureFlow(Dune::Fem::Parameter::getValue<bool>("UseMeanCurvatureFlow",0));
    if(useMeanCurvatureFlow)
      std::cout<<"Problem type: mean curvature flow.\n";
    else
      std::cout<<"Problem type: surface diffusion.\n";

    // compute solution
    typedef Dune::Fem::FemSchemeInterface<GridType> FemSchemeType;
    FemSchemeType femScheme(grid,useMeanCurvatureFlow);
    Dune::Fem::computeInterface<FemSchemeType>(femScheme);

    // dump final mesh as msh
    std::string fileNameFinalMesh(Dune::Fem::Parameter::getValue<std::string>("FileNameFinalMesh",""));
    if(!fileNameFinalMesh.empty())
    {
      const std::string& path(Dune::Fem::Parameter::getValue<std::string>("fem.prefix","."));
      if(!Dune::Fem::directoryExists(path))
        Dune::Fem::createDirectory(path);
      Dune::GmshWriter<typename GridType::LeafGridView> gmshWriter(grid.leafGridView());
      gmshWriter.setPrecision(15);
      gmshWriter.write(path+"/"+fileNameFinalMesh,elementsIDs);
      std::cout<<"\nFinal mesh dumped into "<<fileNameFinalMesh<<".\n";
    }

    // output total running time
    timer.stop();
    std::cout<<"\nTotal running time: "<<timer.elapsed()<<" seconds.\n";

    // output parameters
    std::cout<<"\nParameters used:\n";
    Dune::Fem::Parameter::write(std::cout);

    return 0;
  }

  catch(std::exception& e)
  {
    throw;
  }

  catch(...)
  {
    std::cerr<<"Unknown exception thrown!\n";
    exit(1);
  }
}
