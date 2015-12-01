#ifndef DUNE_FEM_INTERFACEOPERATOR_HH
#define DUNE_FEM_INTERFACEOPERATOR_HH

#include <dune/geometry/referenceelements.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/quadrature/lumpingquadrature.hh>
#include <dune/fem/solver/timeprovider.hh>

#include "normal.hh"

#include <fstream>
#include <vector>
#include <string>

namespace Dune
{
namespace Fem
{

template<class LinearOperatorImp,class TimeProviderImp=FixedStepTimeProvider<>>
class InterfaceOperator:public Operator<typename LinearOperatorImp::DomainFunctionType,typename LinearOperatorImp::RangeFunctionType>
{
  public:
  typedef LinearOperatorImp LinearOperatorType;
  typedef TimeProviderImp TimeProviderType;
  typedef typename LinearOperatorType::DomainFunctionType DomainFunctionType;
  typedef typename LinearOperatorType::RangeFunctionType RangeFunctionType;
  typedef DomainFunctionType DiscreteFunctionType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteSpaceType;
  typedef typename LinearOperatorType::MatrixType MatrixType;
  typedef InterfaceOperator<LinearOperatorType,TimeProviderType> ThisType;

  explicit InterfaceOperator(const DiscreteSpaceType& space,const TimeProviderType& timeProvider,const bool& useMeanCurvFlow):
    space_(space),timeprovider_(timeProvider),op_("interface operator",space_,space_),usemeancurvflow_(useMeanCurvFlow)
  {}

  InterfaceOperator(const ThisType& )=delete;

  virtual inline void operator()(const DomainFunctionType& u,RangeFunctionType& w) const
  {
    op_.apply(u,w);
  }

  // dump system matrix into file
  void print(const std::string& filename="interface_matrix.dat") const
  {
    std::ofstream ofs(filename);
    const auto rows(op_.matrix().rows());
    auto count(decltype(rows){0});
    for(auto row=decltype(rows){0};row!=rows;++row)
    {
      while(count<(op_.matrix().numNonZeros()*(row+1)))
      {
        const auto entry(op_.matrix().realValue(count));
        const auto value(entry.first);
        const auto col(entry.second);
        if((std::abs(value)>1.e-13)&&(col>-1))
          ofs<<row+1<<" "<<col+1<<" "<<value<<std::endl;
        ++count;
      }
    }
  }

  inline const DiscreteSpaceType& domainSpace() const
  {
    return space_;
  }

  inline const DiscreteSpaceType& rangeSpace() const
  {
    return space_;
  }

  inline const LinearOperatorType& systemMatrix() const
  {
    return op_;
  }

  void assemble() const
  {
    // allocate matrix
    DiagonalAndNeighborStencil<DiscreteSpaceType,DiscreteSpaceType> stencil(space_,space_);
    op_.reserve(stencil);
    op_.clear();
    // allocate local basis
    const auto blockSize(DiscreteSpaceType::localBlockSize);
    typedef typename DiscreteFunctionType::LocalFunctionType::RangeType LocalFunctionRangeType;
    std::vector<LocalFunctionRangeType> phi(space_.blockMapper().maxNumDofs()*blockSize );
    typedef typename DiscreteFunctionType::LocalFunctionType::JacobianRangeType LocalFunctionJacobianRangeType;
    std::vector<LocalFunctionJacobianRangeType> gradphi(space_.blockMapper().maxNumDofs()*blockSize);
    // extract dimensions
    constexpr auto worlddim(DiscreteSpaceType::GridType::dimensionworld);
    constexpr auto rangedim(DiscreteSpaceType::FunctionSpaceType::dimRange);
    // define normal
    typedef typename DiscreteSpaceType::RangeFieldType RangeFieldType;
    typedef typename DiscreteSpaceType::GridType::ctype ctype;
    typedef Normal<ctype,worlddim> NormalType;
    NormalType normal;
    typename NormalType::NormalVectorType normalVector;
    // assemble global matrix
    for(const auto entity:space_)
    {
      // compute normal
      normal(entity,normalVector);
      // extract local matrix and basis functions
      auto localMatrix(op_.localMatrix(entity,entity));
      const auto& baseSet(localMatrix.domainBasisFunctionSet());
      // assemble local A_m (curvature) and local \vec{A_m} (position)
      CachingQuadrature<typename DiscreteSpaceType::GridPartType,0> quadrature(entity,2*space_.order()+1);
      for(const auto qp:quadrature)
      {
        // evaluate basis functions and weight
        baseSet.evaluateAll(qp,phi);
        baseSet.jacobianAll(qp,gradphi);
        const auto weight(entity.geometry().integrationElement(qp.position())*qp.weight());
        // fill A_m
        for(auto i=0;i!=worlddim;++i)
        {
          for(auto j=0;j!=worlddim;++j)
          {
            RangeFieldType value(0.0);
            if(usemeancurvflow_)
              value=phi[i][0]*phi[j][0];
            else
              value=gradphi[i][0]*gradphi[j][0];
            value*=weight*timeprovider_.deltaT(); // value=0.0 means v=0
            localMatrix.add(i,j,value);
          }
        }
        // fill \vec{A_m}
        const auto columnLocalSize(localMatrix.columns());
        const auto rowLocalSize(localMatrix.rows());
        for(std::size_t i=worlddim;i!=rowLocalSize;++i)
        {
          for(std::size_t j=worlddim;j!=columnLocalSize;++j)
          {
            RangeFieldType value(0.0);
            for(auto k=1;k!=rangedim;++k)
              value+=gradphi[i][k]*gradphi[j][k];
            value*=weight;
            localMatrix.add(i,j,value);
          }
        }
      }
      // assemble local \vec{N_m} (curvature_j-position_i)
      CachingLumpingQuadrature<typename DiscreteSpaceType::GridPartType,0> lumpingQuadrature(entity);
      for(const auto qp:lumpingQuadrature)
      {
        // evaluate basis functions and weight
        baseSet.evaluateAll(qp,phi);
        const auto weight(entity.geometry().integrationElement(qp.position())*qp.weight());
        // fill \vec{N_m}
        const auto rowLocalSize(localMatrix.rows());
        for(std::size_t i=worlddim;i!=rowLocalSize;++i)
        {
          for(auto j=0;j!=worlddim;++j)
          {
            RangeFieldType value(0.0);
            for(auto index=0;index!=worlddim;++index)
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
  const TimeProviderType& timeprovider_;
  mutable LinearOperatorType op_;
  const bool& usemeancurvflow_;
};

}
}
#endif // DUNE_FEM_INTERFACEOPERATOR_HH
