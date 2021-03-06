#ifndef DUNE_FEM_INTERFACEOPERATOR_HH
#define DUNE_FEM_INTERFACEOPERATOR_HH

#include <dune/fem/io/io.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/quadrature/lumpingquadrature.hh>

#include "normal.hh"

#include <fstream>
#include <string>
#include <vector>

namespace Dune
{
namespace Fem
{

template<typename DiscreteFunctionImp,template<typename ,typename > typename LinearOperatorImp=SparseRowLinearOperator>
class InterfaceOperator:public Operator<DiscreteFunctionImp,DiscreteFunctionImp>
{
  public:
  typedef DiscreteFunctionImp DiscreteFunctionType;
  typedef DiscreteFunctionType DomainFunctionType;
  typedef DiscreteFunctionType RangeFunctionType;
  typedef LinearOperatorImp<DomainFunctionType,RangeFunctionType> LinearOperatorType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteSpaceType;
  typedef typename LinearOperatorType::MatrixType MatrixType;
  typedef InterfaceOperator<DiscreteFunctionType,LinearOperatorImp> ThisType;

  explicit InterfaceOperator(const DiscreteSpaceType& space,bool useMeanCurvFlow):
    space_(space),op_("interface operator",space_,space_),usemeancurvflow_(useMeanCurvFlow)
  {}

  InterfaceOperator(const ThisType& )=delete;

  virtual void operator()(const DomainFunctionType& u,RangeFunctionType& w) const
  {
    op_.apply(u,w);
  }

  void print(const std::string& filename="interface_matrix.dat",unsigned int offset=0) const
  {
    const std::string& path(Parameter::getValue<std::string>("fem.prefix","."));
    if(!directoryExists(path))
      createDirectory(path);
    std::ofstream ofs(path+"/"+filename);
    op_.matrix().print(ofs,offset);
  }

  const DiscreteSpaceType& domainSpace() const
  {
    return space_;
  }

  const DiscreteSpaceType& rangeSpace() const
  {
    return space_;
  }

  const LinearOperatorType& systemMatrix() const
  {
    return op_;
  }

  // assemble operator, use null velocity to compute initial curvature of interface
  template<typename TimeProviderType>
  void assemble(const TimeProviderType& timeProvider,bool velocityNotNull)
  {
    // allocate matrix
    DiagonalAndNeighborStencil<DiscreteSpaceType,DiscreteSpaceType> stencil(space_,space_);
    op_.reserve(stencil);
    op_.clear();
    // allocate local basis
    std::vector<typename DiscreteFunctionType::RangeType> phi(space_.maxNumDofs());
    std::vector<typename DiscreteFunctionType::JacobianRangeType> gradphi(space_.maxNumDofs());
    // extract dimensions
    constexpr unsigned int worlddim(DiscreteSpaceType::GridType::dimensionworld);
    constexpr unsigned int rangedim(DiscreteSpaceType::FunctionSpaceType::dimRange);
    typedef typename DiscreteSpaceType::RangeFieldType RangeFieldType;
    // assemble global matrix
    for(const auto& entity:space_)
    {
      // compute normal
      const auto normalVector(computeNormal(entity));
      // extract local matrix and basis functions
      auto localMatrix(op_.localMatrix(entity,entity));
      const auto columnLocalSize(localMatrix.columns());
      const auto rowLocalSize(localMatrix.rows());
      const auto& baseSet(localMatrix.domainBasisFunctionSet());
      // assemble local matrix
      const CachingLumpingQuadrature<typename DiscreteSpaceType::GridPartType,0> quadrature(entity, 0);
      for(const auto& qp:quadrature)
      {
        // evaluate basis functions and weight
        baseSet.evaluateAll(qp,phi);
        baseSet.jacobianAll(qp,gradphi);
        const auto weight(entity.geometry().integrationElement(qp.position())*qp.weight());
        // fill A_m (curvature)
        if(velocityNotNull)
        {
          for(auto i=decltype(worlddim){0};i!=worlddim;++i)
            for(auto j=decltype(worlddim){0};j!=worlddim;++j)
            {
              RangeFieldType value(0.0);
              if(usemeancurvflow_)
                value=phi[i][0]*phi[j][0];
              else
                value=gradphi[i][0]*gradphi[j][0];
              value*=weight;
              localMatrix.add(i,j,value);
            }
        }
        // fill \vec{A_m} (position)
        for(auto i=decltype(rowLocalSize){worlddim};i!=rowLocalSize;++i)
          for(auto j=decltype(columnLocalSize){worlddim};j!=columnLocalSize;++j)
          {
            RangeFieldType value(0.0);
            for(auto k=decltype(rangedim){1};k!=rangedim;++k)
              value+=gradphi[i][k]*gradphi[j][k];
            value*=weight;
            localMatrix.add(i,j,value);
          }
        // fill \vec{N_m} (curvature_j-position_i) and \vec{N_m}^T
        for(auto i=decltype(rowLocalSize){worlddim};i!=rowLocalSize;++i)
          for(auto j=decltype(worlddim){0};j!=worlddim;++j)
          {
            RangeFieldType value(0.0);
            for(auto index=decltype(worlddim){0};index!=worlddim;++index)
              value+=phi[i][index+1]*normalVector[index];
            value*=weight*phi[j][0];
            localMatrix.add(i,j,value);
            localMatrix.add(j,i,-1.0*value/(timeProvider.deltaT()));
          }
      }
    }
  }

  private:
  const DiscreteSpaceType& space_;
  LinearOperatorType op_;
  const bool usemeancurvflow_;
};

}
}
#endif // DUNE_FEM_INTERFACEOPERATOR_HH
