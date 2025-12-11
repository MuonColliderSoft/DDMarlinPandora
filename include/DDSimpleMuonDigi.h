#ifndef DDSimpleMuonDigi_H
#define DDSimpleMuonDigi_H 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include <vector>

#include "CalorimeterHitType.h"
#include "TFormula.h"
#include "marlin/ProcessorEventSeeder.h"
#include "marlin/AIDAProcessor.h"
#include "marlin/Global.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"


using namespace lcio ;
using namespace marlin ;

namespace EVENT {
  class SimCalorimeterHit ;
}


/** === DDSimpleMuonDigi Processor === <br>
 *  Simple calorimeter digitizer for the muon detectors.
 *  Converts SimCalorimeterHit collections to one 
 *  CalorimeterHit collection applying a threshold and an calibration constant...
 * 
 *  @version $Id$
 */
class DDSimpleMuonDigi : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new DDSimpleMuonDigi ; }
  
  
  DDSimpleMuonDigi() ;
  
  virtual void init() ;
  
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  virtual void end() ;

  bool useLayer(CHT::Layout caloLayout, unsigned int layer) ;
  float computeHitTime( const EVENT::SimCalorimeterHit *h ) const ;
  void smearPosition(const float* pos, float* corr_pos);
  bool timeHitCut(const float* pos, float t);

 protected:

  gsl_rng* _rng {nullptr};

  int _nRun = 0;
  int _nEvt = 0;

  IntVec _layersToKeepBarrelVec{}, _layersToKeepEndcapVec{};
  std::vector<bool>  _useLayersBarrelVec{}, _useLayersEndcapVec{};

  std::vector<std::string> _muonCollections{};

  std::string _outputMuonCollection = "";
  std::string _outputRelCollection = "";

  std::string _cellIDLayerString = "layer";

  float _thresholdMuon = 0.025;
  float _timeThresholdMuon = _thresholdMuon ;
  float _calibrCoeffMuon = 120000;
  float _maxHitEnergyMuon = 2.0;
  
  float _timeResolution = 0.1; // [ns]
  float _cmm = 299.792458;  // [mm/ns]

  std::string _detectorNameBarrel = "YokeBarrel";
  std::string _detectorNameEndcap = "YokeEndcap";
  unsigned int _muonDetBarrel;
  unsigned int _muonDetEndcap;
  
  std::vector<float> _angleRegions{};
  std::vector<float> _timeMins = {-1.0};
  std::vector<float> _timeMaxs = {5.0};
  
} ;

#endif



