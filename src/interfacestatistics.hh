#ifndef DUNE_FEM_MISCDEBUG_HH
#define DUNE_FEM_MISCDEBUG_HH

#include <algorithm>
#include <cstdlib>
#include <limits>

#include "gnuplotwriter.hh"

namespace Dune
{
namespace Fem
{

// dump interface volume
struct InterfaceVolumeInfo:public GnuplotWriter
{
  typedef GnuplotWriter BaseType;

  InterfaceVolumeInfo(unsigned int precision=6):
    BaseType("interface_volume",precision)
  {}

  using BaseType::add;

  template<typename GridPartType,typename TimeProviderType>
  void add(const GridPartType& gridPart,const TimeProviderType& timeProvider)
  {
    double volume(0);
    for(const auto& entity:elements(gridPart))
      volume+=std::abs(entity.geometry().volume());
    add(timeProvider.time(),volume);
  }
};

// dump entity ratio
struct EntityRatioInfo:public GnuplotWriter
{
  typedef GnuplotWriter BaseType;

  EntityRatioInfo(unsigned int precision=6):
    BaseType("entity_ratio",precision)
  {}

  using BaseType::add;

  template<typename GridPartType,typename TimeProviderType>
  void add(const GridPartType& gridPart,const TimeProviderType& timeProvider)
  {
    double minVolume(std::numeric_limits<double>::max());
    double maxVolume(std::numeric_limits<double>::min());
    for(const auto& entity:elements(gridPart))
    {
      const auto volume(std::abs(entity.geometry().volume()));
      minVolume=std::min(volume,minVolume);
      maxVolume=std::max(volume,maxVolume);
    }
    add(timeProvider.time(),maxVolume/minVolume);
  }
};

// dump interface average radius
struct AverageRadiusInfo:public GnuplotWriter
{
  typedef GnuplotWriter BaseType;

  AverageRadiusInfo(unsigned int precision=6):
    BaseType("average_radius",precision)
  {}

  using BaseType::add;

  template<typename GridPartType,typename TimeProviderType>
  void add(const GridPartType& gridPart,const TimeProviderType& timeProvider,
           const typename GridPartType::GridType::template Codim<0>::Entity::Geometry::GlobalCoordinate& center=
           typename GridPartType::GridType::template Codim<0>::Entity::Geometry::GlobalCoordinate(0))
  {
    double radius(0);
    for(const auto& vertex:vertices(gridPart))
    {
      const auto position(vertex.geometry().center()-center);
      radius+=position.two_norm();
    }
    radius/=static_cast<double>(gridPart.grid().size(GridPartType::dimension));
    add(timeProvider.time(),radius);
  }
};

}
}

#endif // DUNE_FEM_MISCDEBUG_HH
