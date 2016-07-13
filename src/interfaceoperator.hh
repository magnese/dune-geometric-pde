#ifndef DUNE_FEM_INTERFACEOPERATOR_HH
#define DUNE_FEM_INTERFACEOPERATOR_HH

#include <dune/geometry/referenceelements.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/quadrature/lumpingquadrature.hh>

#include "normal.hh"

#include <fstream>
#include <vector>
#include <string>

namespace Dune
{
namespace Fem
{

template<typename DiscreteFunctionImp,template<typename ,typename > typename LinearOperatorImp=SparseRowLinearOperator>
class InterfaceOperator:public Operator<DiscreteFunctionImp,DiscreteFunctionImp>
{
  public:
  typedef DiscreteFunctionImp DiscreteFunctionType;
  typedef LinearOperatorImp<DiscreteFunctionType,DiscreteFunctionType> LinearOperatorType;
  typedef DiscreteFunctionType DomainFunctionType;
  typedef DomainFunctionType RangeFunctionType;
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

  void print(const std::string& filename="interface_matrix.dat") const
  {
    std::ofstream ofs(filename);
    op_.matrix().print(ofs);
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

  template<typename TimeProviderType>
  void assemble(const TimeProviderType& timeProvider)
  {
    // allocate matrix
    DiagonalAndNeighborStencil<DiscreteSpaceType,DiscreteSpaceType> stencil(space_,space_);
    op_.reserve(stencil);
    op_.clear();
    // allocate local basis
    constexpr std::size_t blockSize(DiscreteSpaceType::localBlockSize);
    typedef typename DiscreteFunctionType::LocalFunctionType::RangeType LocalFunctionRangeType;
    std::vector<LocalFunctionRangeType> phi(space_.blockMapper().maxNumDofs()*blockSize );
    typedef typename DiscreteFunctionType::LocalFunctionType::JacobianRangeType LocalFunctionJacobianRangeType;
    std::vector<LocalFunctionJacobianRangeType> gradphi(space_.blockMapper().maxNumDofs()*blockSize);
    // extract dimensions
    constexpr unsigned int worlddim(DiscreteSpaceType::GridType::dimensionworld);
    constexpr unsigned int rangedim(DiscreteSpaceType::FunctionSpaceType::dimRange);
    // define normal
    typedef typename DiscreteSpaceType::RangeFieldType RangeFieldType;
    typedef typename DiscreteSpaceType::GridType::ctype ctype;
    typedef Normal<ctype,worlddim> NormalType;
    NormalType normal;
    typename NormalType::NormalVectorType normalVector;
    // assemble global matrix
    for(const auto& entity:space_)
    {
      // compute normal
      normal(entity,normalVector);
      // extract local matrix and basis functions
      auto localMatrix(op_.localMatrix(entity,entity));
      const auto& baseSet(localMatrix.domainBasisFunctionSet());
      // assemble local A_m (curvature) and local \vec{A_m} (position)
      CachingQuadrature<typename DiscreteSpaceType::GridPartType,0> quadrature(entity,2*space_.order()+1);
      for(const auto& qp:quadrature)
      {
        // evaluate basis functions and weight
        baseSet.evaluateAll(qp,phi);
        baseSet.jacobianAll(qp,gradphi);
        const auto weight(entity.geometry().integrationElement(qp.position())*qp.weight());
        // fill A_m
        for(auto i=decltype(worlddim){0};i!=worlddim;++i)
        {
          for(auto j=decltype(worlddim){0};j!=worlddim;++j)
          {
            RangeFieldType value(0.0);
            if(usemeancurvflow_)
              value=phi[i][0]*phi[j][0];
            else
              value=gradphi[i][0]*gradphi[j][0];
            value*=weight*timeProvider.deltaT(); // value=0.0 means v=0
            localMatrix.add(i,j,value);
          }
        }
        // fill \vec{A_m}
        const auto columnLocalSize(localMatrix.columns());
        const auto rowLocalSize(localMatrix.rows());
        for(auto i=decltype(rowLocalSize){worlddim};i!=rowLocalSize;++i)
        {
          for(auto j=decltype(columnLocalSize){worlddim};j!=columnLocalSize;++j)
          {
            RangeFieldType value(0.0);
            for(auto k=decltype(rangedim){1};k!=rangedim;++k)
              value+=gradphi[i][k]*gradphi[j][k];
            value*=weight;
            localMatrix.add(i,j,value);
          }
        }
      }
      // assemble local \vec{N_m} (curvature_j-position_i)
      CachingLumpingQuadrature<typename DiscreteSpaceType::GridPartType,0> lumpingQuadrature(entity);
      for(const auto& qp:lumpingQuadrature)
      {
        // evaluate basis functions and weight
        baseSet.evaluateAll(qp,phi);
        const auto weight(entity.geometry().integrationElement(qp.position())*qp.weight());
        // fill \vec{N_m}
        const auto rowLocalSize(localMatrix.rows());
        for(auto i=decltype(rowLocalSize){worlddim};i!=rowLocalSize;++i)
        {
          for(auto j=decltype(worlddim){0};j!=worlddim;++j)
          {
            RangeFieldType value(0.0);
            for(auto index=decltype(worlddim){0};index!=worlddim;++index)
              value+=phi[i][index+1]*normalVector[index];
            value*=weight*phi[j][0];
            localMatrix.add(i,j,value);
            localMatrix.add(j,i,-1.0*value);
          }
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
