#ifndef DUNE_FEM_ASSEMBLEINTERFACERHS_HH
#define DUNE_FEM_ASSEMBLEINTERFACERHS_HH

namespace Dune
{
namespace Fem
{

template<typename DiscreteFunctionType,typename OperatorType>
void assembleInterfaceRHS(DiscreteFunctionType& rhs,const OperatorType& op)
{
  rhs.clear();

  DiscreteFunctionType temp("temp",rhs.space());
  temp.clear();
  temp.template subDiscreteFunction<1>()-=rhs.space().grid().coordFunction().discreteFunction();

  op(temp,rhs);

  rhs.template subDiscreteFunction<0>().clear();
}

}
}

#endif // DUNE_FEM_ASSEMBLEINTERFACERHS_HH
