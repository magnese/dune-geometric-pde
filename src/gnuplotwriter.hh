#ifndef DUNE_FEM_GNUPLOTWRITER_HH
#define DUNE_FEM_GNUPLOTWRITER_HH

#include <fstream>
#include <iomanip>
#include <list>
#include <string>
#include <tuple>
#include <utility>

#include <dune/fem/io/parameter.hh>

namespace Dune
{
namespace Fem
{

// generic gnuplot writer
struct GnuplotWriter
{
  typedef std::list<std::tuple<double,double>> ListType;

  GnuplotWriter(const std::string& fileName,unsigned int precision=6):
    filename_(Parameter::getValue<std::string>("fem.prefix",".")+"/"+fileName+".dat"),precision_(precision)
  {}

  void add(double first,double second)
  {
    values_.emplace_back(first,second);
  }

  bool isEmpty() const
  {
    return values_.size()==0;
  }

  ~GnuplotWriter()
  {
    finalize();
  }

  void finalize() const
  {
    if(!isEmpty())
    {
      std::ofstream file(filename_);
      file<<std::setprecision(precision_);
      for(const auto& value:values_)
        file<<std::get<0>(value)<<" "<<std::get<1>(value)<<std::endl;
    }
  }

  std::string filename_;
  unsigned int precision_;
  ListType values_;
};

}
}

#endif // DUNE_FEM_GNUPLOTWRITER_HH
