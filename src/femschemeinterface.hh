#ifndef DUEN_FEM_FEMSCHEMEINTERFACE_HH
#define DUEN_FEM_FEMSCHEMEINTERFACE_HH

#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/solver/umfpacksolver.hh>
#include <dune/fem/solver/timeprovider.hh>

#include "interfaceoperator.hh"
#include "interfacerhs.hh"

namespace Dune
{
namespace Fem
{

template<class GridImp,class TimeProviderImp=FixedStepTimeProvider<>>
class FemSchemeInterface
{
  public:
  // define grid types
  typedef GridImp GridType;
  typedef FemSchemeInterface<GridType> ThisType;
  typedef LeafGridPart<GridType> GridPartType;

  // define spaces and functions
  typedef FunctionSpace<double,double,GridType::dimensionworld,1> CurvatureContinuosSpaceType;
  typedef FunctionSpace<double,double,GridType::dimensionworld,GridType::dimensionworld> PositionContinuosSpaceType;
  typedef LagrangeDiscreteFunctionSpace<CurvatureContinuosSpaceType,GridPartType,POLORDER> CurvatureDiscreteSpaceType;
  typedef LagrangeDiscreteFunctionSpace<PositionContinuosSpaceType,GridPartType,POLORDER> PositionDiscreteSpaceType;
  typedef TupleDiscreteFunctionSpace<CurvatureDiscreteSpaceType,PositionDiscreteSpaceType> CombinedDiscreteSpaceType;
  typedef AdaptiveDiscreteFunction<CombinedDiscreteSpaceType> CombinedDiscreteFunctionType;

  // define operator and rhs
  typedef SparseRowLinearOperator<CombinedDiscreteFunctionType,CombinedDiscreteFunctionType> LinearOperatorType;
  typedef InterfaceOperator<LinearOperatorType> InterfaceOperatorType;
  typedef InterfaceRHS<CombinedDiscreteFunctionType> InterfaceRHSType;

  // define inverse operator
  typedef UMFPACKOp<CombinedDiscreteFunctionType,InterfaceOperatorType> InterfaceInverseOperatorType;

  // define time provider
  typedef TimeProviderImp TimeProviderType;

  explicit FemSchemeInterface(GridType& grid,const bool& useMeanCurvFlow):
    grid_(grid),gridpart_(grid_),space_(gridpart_),usemeancurvflow_(useMeanCurvFlow)
  {}

  FemSchemeInterface(const ThisType& )=delete;

  inline GridType& getGrid()
  {
    return grid_;
  }
  inline const GridPartType& getGridPart() const
  {
    return gridpart_;
  }
  inline const CombinedDiscreteSpaceType& getSpace() const
  {
    return space_;
  }

  // setup and solve the linear system
  void operator()(CombinedDiscreteFunctionType& solution,const TimeProviderType& timeProvider)
  {
    // clear solution
    solution.clear();
    // assemble operator
    InterfaceOperatorType op(space_,timeProvider,usemeancurvflow_);
    op.assemble();
    // assemble rhs
    InterfaceRHSType RHS(space_);
    RHS.assemble(op);
    // solve the linear system
    InterfaceInverseOperatorType interfaceInvOp(op);
    interfaceInvOp(RHS.rhs(),solution);
  }

  private:
  GridType& grid_;
  GridPartType gridpart_;
  CombinedDiscreteSpaceType space_;
  const bool& usemeancurvflow_;
};

}
}

#endif // DUEN_FEM_FEMSCHEMEINTERFACE_HH
