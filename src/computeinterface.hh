#ifndef DUNE_FEM_COMPUTEINTERFACE_HH
#define DUNE_FEM_COMPUTEINTERFACE_HH

#include <iostream>
#include <tuple>
#include <cmath>

#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/io/file/datawriter.hh>
#include <dune/common/timer.hh>
#include <dune/fem/io/parameter.hh>

namespace Dune
{
namespace Fem
{

template<class FemSchemeType>
void computeInterface(FemSchemeType& femScheme)
{
  // create time provider
  FixedStepTimeProvider<> timeProvider;

  // get grid
  auto& grid(femScheme.getGrid());

  // get space
  typedef typename FemSchemeType::CombinedDiscreteFunctionType CombinedDiscreteFunctionType;
  CombinedDiscreteFunctionType solution("solution",femScheme.getSpace());
  solution.clear();

  // output number of dofs
  const auto numDofsCurvature(femScheme.getSpace().template subDiscreteFunctionSpace<0>().size());
  const auto numDofsPosition(femScheme.getSpace().template subDiscreteFunctionSpace<1>().size());
  std::cout<<"Solving for "<<numDofsCurvature<<" unkowns for curvature and "<<numDofsPosition<<" unkowns for position."<<std::endl;

  // create structure to dump on file
  typedef std::tuple<const CombinedDiscreteFunctionType*> IOTupleType;
  IOTupleType ioTuple(&solution);
  typedef DataOutput<typename FemSchemeType::GridType,IOTupleType> DataOutputType;
  DataOutputType dataOutput(grid,ioTuple);

  // dump bulk solution at t0 and advance time provider
  dataOutput.write(timeProvider);
  timeProvider.next();

  // extract iterators which point to the first and to the last dof of position
  const auto positionItBegin(std::next(solution.dbegin(),numDofsCurvature));
  const auto positionItEnd(solution.dend());

  // enable/disable check interface is stationary
  bool interfaceStationary(true);
  const bool createStationaryInterface(Parameter::getValue<bool>("CreateStationaryInterface",0));
  if(createStationaryInterface)
    std::cout<<std::endl<<"WARNING: the scheme will run until the interface is stationary!"<<std::endl;

  // solve
  const double endTime(Parameter::getValue<double>("EndTime",1.0)+0.1*timeProvider.deltaT());
  for(;(timeProvider.time()<=endTime)||(!interfaceStationary);timeProvider.next())
  {
    // print time
    std::cout<<std::endl<<"Time step "<<timeProvider.timeStep()<<" (time = "<<timeProvider.time()<<" s)."<<std::endl;
    // start timer
    Timer timer(false);
    timer.start();
    // compute solution
    femScheme(solution,timeProvider);
    // check if the interface is stationary
    if(createStationaryInterface)
    {
      interfaceStationary=true;
      for(auto positionIt=positionItBegin;positionIt!=positionItEnd&&interfaceStationary;++positionIt)
        if(std::abs(*positionIt)>1.e-15)
          interfaceStationary=false;
      if(interfaceStationary)
        std::cout<<"Interface is stationary."<<std::endl;
      else
        std::cout<<"Interface is NOT stationary."<<std::endl;
    }
    // update grid and solution
    auto coordIt(grid.coordFunction().discreteFunction().dbegin());
    for(auto positionIt=positionItBegin;positionIt!=positionItEnd;++positionIt,++coordIt)
    {
      (*coordIt)+=(*positionIt); // X+dX
      (*positionIt)=(*coordIt); // dump on file X instead of dX
    }
    // stop timer
    timer.stop();
    std::cout<<"Time elapsed for assembling and solving : "<<timer.elapsed()<<" seconds."<<std::endl;
    // dump solution on file
    dataOutput.write(timeProvider);
  }
}

}
}

#endif // DUNE_FEM_COMPUTEINTERFACE_HH
