//  Telescoping Jets Package
//
//  Alexander Emerman, Yang-Ting Chien, Shih-Chieh Hsu
//

#include "TelescopingJets.hh"

  // Returns the mass volatility for N axes
double TelescopingJets::result(int N, const fastjet::PseudoJet& jet) const {

  int nconst = jet.constituents().size();
  //std::cout<<"Number of constits: "<<nconst<<std::endl;
  if(nconst==1){
    //std::cout<<"Only one constituent, can't calculate TJets("<<N<<")"<<std::endl;
    return -1;
  }

  std::vector<double> tjet_masses = getTjetMasses(N,jet);
  return getVolatility(tjet_masses);
}
  // Returns N subjets per R value
std::vector<std::vector<fastjet::PseudoJet> > TelescopingJets::getSubjets(int N, const fastjet::PseudoJet& jet) const {
  _setSubjets(N, jet);
  return _current_subjets;
}
  // Returns the sum of the subjets at each R value
std::vector<fastjet::PseudoJet> TelescopingJets::getTjets(int N, const fastjet::PseudoJet& jet) const {
  _setTjets(N, jet);
  return _current_tjets;
}
  // Return the invariant mass of the subjets at each R value
std::vector<double> TelescopingJets::getTjetMasses(int N, const fastjet::PseudoJet& jet) const {
  _setTjetMasses(N, jet);
  return _current_tjet_masses;
}

 // Useful mathematical functions to find the average, rms and volatility
 // of a set of values
double getAverage(const std::vector<double>& values) {
  if (values.size() == 0) throw std::length_error("Asking for average of empty vector.\n");
  double v = 0;
  for (unsigned int i=0; i < values.size(); i++) {
    v += values.at(i);
  }
  return v/values.size();
}
double getRMS(const std::vector<double>& values) {
  if (values.size() == 0) throw std::length_error("Asking for rms of empty vector.\n");
  double v = getAverage(values);
  double v2 = 0;
  for(unsigned int i=0; i < values.size(); i++) {
      v2 += (values.at(i) - v) * (values.at(i) - v);
  }
  double sqRMS = v2/values.size();
  return std::sqrt(sqRMS);
}
double getVolatility(const std::vector<double>& values) {
  if (values.size() == 0) throw std::length_error("Asking for volatility of empty vector.\n");
  return getRMS(values)/getAverage(values);
}

  // set user-defined axes
  // currently does nothing because axes are reset before every calculation
// void TelescopingJets::setAxes(const std::vector<fastjet::PseudoJet>& my_axes) {
//   _N = my_axes.size();
//   _current_axes = my_axes;
//   _axes_are_set = true;
//   return;
// }

//**************************************
// Private methods to do the actual work

void TelescopingJets::_setSubjets(int N, const fastjet::PseudoJet& jet) const {
  _setAxes(N, jet.constituents());
  _setSubjets(jet);
  return;
}
void TelescopingJets::_setSubjets(const fastjet::PseudoJet& jet) const {
  if (!_axes_are_set) throw std::runtime_error("axes not set");
  // reset objects so that they can never hold information on a set of particles other than the current one
  _current_tjet_masses.clear();
  _current_tjets.clear();
  _current_subjets.clear();
  _current_rejected.clear();

  std::vector<fastjet::PseudoJet> particles = jet.constituents();
  for (unsigned int i=0; i < _r_values.size(); i++) {
    // partition particles by nearest axis
    std::vector<std::vector<fastjet::PseudoJet> > jet_partition(_N+1);
    _partitionJet(particles,_r_values.at(i),jet_partition);

    std::vector<fastjet::PseudoJet> subjet_vec(_N, fastjet::PseudoJet(0.,0.,0.,0.));
    for (unsigned int n=0; n < _N; n++) {
      // use the join method so that subjets have valid constituent structure
      subjet_vec.at(n) = join( jet_partition.at(n) );
    } // end loop over axes

    _current_subjets.push_back(subjet_vec);
    _current_rejected.push_back( jet_partition.at(_N) );
  } // end loop over R values
  return;
}

void TelescopingJets::_setTjets(int N, const fastjet::PseudoJet&  jet) const {
  _setAxes(N, jet.constituents());
  _setTjets(jet);
  return;
}
void TelescopingJets::_setTjets(const fastjet::PseudoJet&  jet) const {
  if (!_axes_are_set) throw std::runtime_error("axes not set");
  _setSubjets(jet);
  // reset objects so that they can never hold information on a set of particles other than the current one
  _current_tjet_masses.clear();
  _current_tjets.clear();

  for (unsigned int i=0; i < _r_values.size(); i++) {
    std::vector<fastjet::PseudoJet> temp_vec;
    for (unsigned int j=0; j < _N; j++) {
      temp_vec.push_back( (_current_subjets.at(i)).at(j) );
    } // end loop over axes
    // use the join method so that t-jets have valid constituent structure
    _current_tjets.push_back(join(temp_vec));
  } // end loop over R values
  return;
}

void TelescopingJets::_setTjetMasses(int N, const fastjet::PseudoJet& jet) const {
  _setAxes(N, jet.constituents());
  _setTjetMasses(jet);
  return;
}
void TelescopingJets::_setTjetMasses(const fastjet::PseudoJet& jet) const {
  if (!_axes_are_set) throw std::runtime_error("axes not set");
  _setTjets(jet);
  // reset objects so that they can never hold information on a set of particles other than the current one
  _current_tjet_masses.clear();
  for (unsigned int i=0; i < _r_values.size(); i++) {
    _current_tjet_masses.push_back( (_current_tjets.at(i)).m() );
  } // end loop over R values
  return;
}

void TelescopingJets::_setAxes(int N, const std::vector<fastjet::PseudoJet>& particles) const {
  _N = N;

/* OLD CODE
  std::vector<fastjet::PseudoJet> seed_axes;
  seed_axes.resize(_N, fastjet::PseudoJet(0.,0.,0.,0.) );
  _current_axes = _axes_finder->getAxes(_N, particles, seed_axes);
  _axes_are_set = true;
*/

  _current_axes = _axes_def->get_axes(_N, particles);
  _axes_are_set = true;

  return;
}

void TelescopingJets::_partitionJet(const std::vector<fastjet::PseudoJet>& particles, double R, std::vector<std::vector<fastjet::PseudoJet> >& jet_partition) const {
  if (_current_axes.size() != _N) throw std::logic_error("Subjet axes uninitialized.\n");
  if (jet_partition.size() != _N+1) throw std::logic_error("Partition uninitialized.\n");

  for (unsigned int i=0; i < particles.size(); i++) {
    double delR, delR_min(999.);
    int j_min(0);
    for (unsigned int j=0; j < _N; j++) {
      delR = (_current_axes.at(j)).delta_R(particles.at(i));
      if (delR < delR_min) {
        delR_min = delR;
        j_min = j;
      }
    } // end loop over axes
    if (delR_min < R) {
      (jet_partition.at(j_min)).push_back(particles.at(i));
    } else {
      (jet_partition.at(_N)).push_back(particles.at(i)); // extra pieces
    }
  }
  return;
}
