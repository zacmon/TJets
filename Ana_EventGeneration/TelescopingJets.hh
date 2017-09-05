//  Telescoping Jets Package
//
//  Alexander Emerman, Yang-Ting Chien, Shih-Chieh Hsu
//

#ifndef __TELESCOPINGJETS_HH__
#define __TELESCOPINGJETS_HH__

//#include "fastjet/contrib/AxesFinder.hh"
#include "fastjet/contrib/AxesDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/FunctionOfPseudoJet.hh"

#include <cmath>
#include <vector>
#include <stdexcept>

// useful math functions in case the user just wants a final
// number rather than dealing with vectors of PseudoJets
double getVolatility(const std::vector<double>& values);
double getAverage(const std::vector<double>& values);
double getRMS(const std::vector<double>& values);

//------------------------------------------------------------------------
/// \class TelescopingJets
//
class TelescopingJets {
public:
  // Main Constructor
  TelescopingJets(const fastjet::contrib::AxesDefinition& axes_def, const std::vector<double> r_values)
  : _r_values(r_values), _axes_def(axes_def.create()) {
    // Don't use AxesDefinitions that require MeasureFunctions
    fastjet::contrib::UnnormalizedMeasure measure(1.); // need a measure function to create the axes finder, it shouldn't be used though
    //_axes_finder.reset(_axes_def->createStartingAxesFinder(measure));
  }
  // Constructor from user axes
//  TelescopingJets(const std::vector<fastjet::PseudoJet> myAxes,
//                  const std::vector<double> r_values)
//  : _r_values(r_values), _axes_def(NULL), _axes_finder(NULL) {
//    setAxes(myAxes);
//  }

  std::string description() const { return "Telescoping jets calculator"; }

  // Methods that do the calculation and return values
    // Returns the mass volatility for N axes
  double result(int N, const fastjet::PseudoJet& jet) const;
  double operator() (int N, const fastjet::PseudoJet& jet) const {
    return result(N,jet);
  }
    // Returns N subjets per R value
    // note: subjets and t-jets should have valid constituent information
  std::vector<std::vector<fastjet::PseudoJet> > getSubjets(int N, const fastjet::PseudoJet& jet) const;
    // Returns the sum of the subjets at each R value
    // TODO: allow different recombination schemes
  std::vector<fastjet::PseudoJet> getTjets(int N, const fastjet::PseudoJet& jet) const;
    // Return the invariant mass of the subjets at each R value
  std::vector<double> getTjetMasses(int N, const fastjet::PseudoJet& jet) const;


  // Methods that return previously calculated values
  std::vector<std::vector<fastjet::PseudoJet> > getCurrentSubjets() const {
    if (_current_subjets.size() == 0) throw std::logic_error("Subjets asked for but not set.\n");
    return _current_subjets;
  }
  std::vector<fastjet::PseudoJet> getCurrentTjets() const {
    if (_current_tjets.size() == 0) throw std::logic_error("Tjets asked for but not set.\n");
    return _current_tjets;
  }
  std::vector<double> getCurrentTjetMasses() const {
    if (_current_tjet_masses.size() == 0) throw std::logic_error("Tjet masses asked for but not set.\n");
    return _current_tjet_masses;
  }
  std::vector<std::vector<fastjet::PseudoJet> > getCurrentRejected() const {
    if (_current_rejected.size() == 0) throw std::logic_error("Rejected pieces asked for but not set.\n");
    return _current_rejected;
  }

  void setAxes(const std::vector<fastjet::PseudoJet>& myAxes);

private:

  std::vector<double> _r_values;
  fastjet::SharedPtr<const fastjet::contrib::AxesDefinition> _axes_def;
//  fastjet::SharedPtr<const fastjet::contrib::AxesFinder> _axes_finder;

  // the calculation methods can only change these values
  mutable unsigned int _N; // current number of axes
  mutable std::vector<fastjet::PseudoJet> _current_axes;
  mutable bool _axes_are_set;
  mutable std::vector<std::vector<fastjet::PseudoJet> > _current_subjets;
  mutable std::vector<std::vector<fastjet::PseudoJet> > _current_rejected;
  mutable std::vector<fastjet::PseudoJet> _current_tjets;
  mutable std::vector<double> _current_tjet_masses;

  // functions that do the actual work. not user accessible
  void _setAxes(int N, const std::vector<fastjet::PseudoJet>& particles) const;

  void _setSubjets(int N, const fastjet::PseudoJet& jet) const;
  void _setSubjets(const fastjet::PseudoJet& jet) const;
  void _setTjets(int N, const fastjet::PseudoJet& jet) const;
  void _setTjets(const fastjet::PseudoJet& jet) const;
  void _setTjetMasses(int N, const fastjet::PseudoJet& jet) const;
  void _setTjetMasses(const fastjet::PseudoJet& jet) const;
  void _partitionJet(const std::vector<fastjet::PseudoJet>& particles, double R, std::vector<std::vector<fastjet::PseudoJet> >& jet_partition) const;

};

#endif  // __TELESCOPINGJETS_HH__
