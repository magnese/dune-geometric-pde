#ifndef DUEN_FEM_FEMSCHEMEINTERFACE_HH
#define DUEN_FEM_FEMSCHEMEINTERFACE_HH

#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/function/tuplediscretefunction.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/solver/umfpacksolver.hh>
#include <dune/fem/solver/spqrsolver.hh>

#include "interfaceoperator.hh"
#include "assembleinterfacerhs.hh"

namespace Dune
{
namespace Fem
{

template<typename GridImp>
class FemSchemeInterface
{
  public:
  // define grid types
  typedef GridImp GridType;
  typedef FemSchemeInterface<GridType> ThisType;
  typedef LeafGridPart<GridType> GridPartType;

  // define spaces and functions
  typedef FunctionSpace<double,double,GridType::dimensionworld,1> CurvatureContinuosSpaceType;
  typedef FunctionSpace<double,double,GridType::dimensionworld,GridType::dimensionworld> DisplacementContinuosSpaceType;
  typedef LagrangeDiscreteFunctionSpace<CurvatureContinuosSpaceType,GridPartType,POLORDER> CurvatureDiscreteSpaceType;
  typedef LagrangeDiscreteFunctionSpace<DisplacementContinuosSpaceType,GridPartType,POLORDER> DisplacementDiscreteSpaceType;
  typedef AdaptiveDiscreteFunction<CurvatureDiscreteSpaceType> CurvatureDiscreteFunctionType;
  typedef AdaptiveDiscreteFunction<DisplacementDiscreteSpaceType> DisplacementDiscreteFunctionType;
  typedef TupleDiscreteFunction<CurvatureDiscreteFunctionType,DisplacementDiscreteFunctionType> DiscreteFunctionType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteSpaceType;

  // define operator
  typedef InterfaceOperator<DiscreteFunctionType> InterfaceOperatorType;

  // define inverse operator
  #if SOLVER_TYPE == 0
  typedef UMFPACKOp<DiscreteFunctionType,InterfaceOperatorType> InterfaceInverseOperatorType;
  #elif SOLVER_TYPE == 1
  typedef SPQROp<DiscreteFunctionType,InterfaceOperatorType> InterfaceInverseOperatorType;
  #endif

  explicit FemSchemeInterface(GridType& grid,bool useMeanCurvFlow):
    grid_(grid),gridpart_(grid_),space_(gridpart_),usemeancurvflow_(useMeanCurvFlow)
  {}

  FemSchemeInterface(const ThisType& )=delete;

  GridType& grid()
  {
    return grid_;
  }
  const GridPartType& gridPart() const
  {
    return gridpart_;
  }
  const DiscreteSpaceType& space() const
  {
    return space_;
  }

  // compute intial curvature
  template<typename TimeProviderType>
  void computeInitialCurvature(DiscreteFunctionType& solution,const TimeProviderType& timeProvider)
  {
    operator()(solution,timeProvider,false);
  }

  // compute solution
  template<typename TimeProviderType>
  void operator()(DiscreteFunctionType& solution,const TimeProviderType& timeProvider,bool velocityNotNull=true)
  {
    // clear solution
    solution.clear();
    // assemble operator
    InterfaceOperatorType op(space_,usemeancurvflow_);
    op.assemble(timeProvider,velocityNotNull);
    // assemble rhs
    DiscreteFunctionType rhs("interface RHS",space_);
    assembleInterfaceRHS(rhs,op);
    // solve the linear system
    InterfaceInverseOperatorType interfaceInvOp(op);
    interfaceInvOp(rhs,solution);
  }

  private:
  GridType& grid_;
  GridPartType gridpart_;
  const DiscreteSpaceType space_;
  const bool usemeancurvflow_;
};

}
}

#endif // DUEN_FEM_FEMSCHEMEINTERFACE_HH
