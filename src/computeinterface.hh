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

template<typename FemSchemeType>
void computeInterface(FemSchemeType& femScheme)
{
  // create time provider
  FixedStepTimeProvider<> timeProvider;

  // get grid
  auto& grid(femScheme.grid());

  // create solution
  typedef typename FemSchemeType::DiscreteFunctionType DiscreteFunctionType;
  DiscreteFunctionType solution("solution",femScheme.space());
  auto& curvature(solution.template subDiscreteFunction<0>());
  curvature.name()="curvature";
  auto& displacement(solution.template subDiscreteFunction<1>());
  displacement.name()="displacement";
  std::cout<<"Solving for "<<curvature.size()<<" unkowns for "<<curvature.name()<<" and "<<displacement.size()<<" unkowns for "<<
    displacement.name()<<"\n";

  // create structure to dump on file
  auto ioTuple(std::make_tuple(&curvature));
  DataOutput<typename FemSchemeType::GridType,decltype(ioTuple)> dataOutput(grid,ioTuple);

  // dump bulk solution at t0 and advance time provider
  femScheme.computeInitialCurvature(solution,timeProvider);
  dataOutput.write(timeProvider);
  timeProvider.next();

  // enable/disable check interface is stationary
  bool interfaceStationary(true);
  const bool createStationaryInterface(Parameter::getValue<bool>("CreateStationaryInterface",0));
  if(createStationaryInterface)
    std::cout<<"\nWARNING: the scheme will run until the interface is stationary!\n";

  // solve
  const double endTime(Parameter::getValue<double>("EndTime",1.0)+0.1*timeProvider.deltaT());
  for(;(timeProvider.time()<=endTime)||(!interfaceStationary);timeProvider.next())
  {
    // print time
    std::cout<<"\nTime step "<<timeProvider.timeStep()<<" (time = "<<timeProvider.time()<<" s).\n";
    // start timer
    Timer timer(false);
    timer.start();
    // compute solution
    femScheme(solution,timeProvider);
    // check if the interface is stationary
    if(createStationaryInterface)
    {
      interfaceStationary=true;
      for(const auto& dof:dofs(displacement))
        if(std::abs(dof)>1.e-15)
        {
          interfaceStationary=false;
          break;
        }
      if(interfaceStationary)
        std::cout<<"Interface is stationary.\n";
      else
        std::cout<<"Interface is NOT stationary.\n";
    }
    // update grid
    grid.coordFunction()+=displacement;
    // stop timer
    timer.stop();
    std::cout<<"Time elapsed for assembling and solving : "<<timer.elapsed()<<" seconds.\n";
    // dump solution on file
    dataOutput.write(timeProvider);
  }
}

}
}

#endif // DUNE_FEM_COMPUTEINTERFACE_HH
