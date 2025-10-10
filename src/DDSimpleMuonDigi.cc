#include "DDSimpleMuonDigi.h"
#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <IMPL/CalorimeterHitImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCFlagImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/LCParameters.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCRelationNavigator.h>

// #include <algorithm>
// #include <string>
#include <cctype> 
#include <cstdlib>  // abs

#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/DetectorData.h"

#include "math.h"

using namespace std;
using namespace lcio ;
using namespace marlin ;


DDSimpleMuonDigi aDDSimpleMuonDigi ;


DDSimpleMuonDigi::DDSimpleMuonDigi() : Processor("DDSimpleMuonDigi") {

  _description = "Performs simple digitization of sim muon hits..." ;
  
  std::vector<std::string> muonCollections;

  muonCollections.push_back(std::string("yoke03_MuonBarrel"));
  muonCollections.push_back(std::string("yoke03_MuonEndCap"));
  muonCollections.push_back(std::string("yoke03_MuonPlug"));
    

  registerInputCollections( LCIO::SIMCALORIMETERHIT, 
			    "MUONCollections" , 
			    "Muon Collection Names" ,
			    _muonCollections ,
			    muonCollections);
  
  registerOutputCollection( LCIO::CALORIMETERHIT, 
			    "MUONOutputCollection" , 
			    "Muon Collection of real Hits" , 
			    _outputMuonCollection , 
			    std::string("MUON")) ; 
  
  registerOutputCollection( LCIO::LCRELATION, 
			    "RelationOutputCollection" , 
			    "CaloHit Relation Collection" , 
			    _outputRelCollection , 
			    std::string("RelationMuonHit")) ; 
  
  registerProcessorParameter("MuonThreshold" , 
			     "Threshold for Muon Hits in GeV" ,
			     _thresholdMuon,
			     (float)0.025);

    registerProcessorParameter("MuonTimeThreshold" ,
			     "Energy threshold for timing information for Muon Hits in GeV" ,
			     _timeThresholdMuon,
			     (float)0.025);

  registerProcessorParameter("CalibrMUON" , 
			     "Calibration coefficients for MUON" ,
			     _calibrCoeffMuon,
			     (float)120000.);

  registerProcessorParameter("MaxHitEnergyMUON", 
			     "maximum hit energy for a MUON hit" ,
			     _maxHitEnergyMuon,
			     (float)2.0);

  IntVec keepBarrelLayersVec, keepEndcapLayersVec;

  registerProcessorParameter("KeepBarrelLayersVec" , 
			     "Vector of Barrel layers to be kept. Layers start at 1!",
			     _layersToKeepBarrelVec,
			     keepBarrelLayersVec);

  registerProcessorParameter("KeepEndcapLayersVec" , 
			     "Vector of Endcap layers to be kept. Layers start at 1!",
			     _layersToKeepEndcapVec,
			     keepEndcapLayersVec);


  registerProcessorParameter("CellIDLayerString" ,
			     "name of the part of the cellID that holds the layer" , 
			     _cellIDLayerString , 
			     std::string("layer")
			     );
  
  registerProcessorParameter("DetectorNameBarrel" ,
                             "Name of the barrel subdetector" , 
                             _detectorNameBarrel , 
                             std::string("YokeBarrel"));
  
  registerProcessorParameter("DetectorNameEndcap" ,
                             "Name of the endcap subdetector" , 
                             _detectorNameEndcap, 
                             std::string("YokeEndcap"));

  //LP: new parameters

  registerProcessorParameter( "MuonHitsTimeResolution",
			      "Time resolution of the muon detector hits",
			      _timeResolution,
			      _timeResolution
			      );

  registerProcessorParameter( "HitAngleRegions",
			      "Hit angle region delimiters",
			      _angleRegions,
			      _angleRegions
			      );

  registerProcessorParameter( "MinTimeHitCuts",
			      "Min time cut by hit region",
			      _timeMins,
			      _timeMins
			      );

  registerProcessorParameter( "MaxTimeHitCuts",
			      "Max time cut by hit region",
			      _timeMaxs,
			      _timeMaxs
			      );
}

void DDSimpleMuonDigi::init() {

  _nRun = -1;
  _nEvt = 0;

  //fg: there cannot be any reasonable default for this string - so we set it to sth. that will cause an exception in case 
  //    the cellID encoding string is not in the collection: 
  UTIL::CellIDDecoder<CalorimeterHit>::setDefaultEncoding("undefined_cellID_encoding:100");


  //Get the number of Layers in the Endcap
  int layersEndcap=0, layersBarrel=0;

  try{
    
    dd4hep::Detector &mainDetector = dd4hep::Detector::getInstance();
    dd4hep::DetElement theDetector = mainDetector.detector(_detectorNameBarrel);
    const dd4hep::rec::LayeredCalorimeterData *yokeBarrelParameters  = theDetector.extension<dd4hep::rec::LayeredCalorimeterData>();
    layersBarrel =  yokeBarrelParameters->layers.size();
  }catch( std::exception& e ){
    streamlog_out(WARNING) << "  oops - no Yoke Barrel available: "<<e.what() << std::endl ;
  }
  try{
    
    dd4hep::Detector &mainDetector = dd4hep::Detector::getInstance();
    dd4hep::DetElement theDetector = mainDetector.detector(_detectorNameEndcap);
    const dd4hep::rec::LayeredCalorimeterData *yokeEndcapParameters  = theDetector.extension<dd4hep::rec::LayeredCalorimeterData>();
    layersEndcap =  yokeEndcapParameters->layers.size();
  }catch(std::exception& e ){
    streamlog_out(WARNING) << "  oops - no Yoke Endcap available: "<<e.what() << std::endl ;
  }

  //If the vectors are empty, we are keeping everything 
  if(_layersToKeepBarrelVec.size() > 0) {
    //layers start at 0
    for(int i = 0; i < layersBarrel; ++i) {
      _useLayersBarrelVec.push_back(false);
      for(IntVec::iterator iter = _layersToKeepBarrelVec.begin(); iter < _layersToKeepBarrelVec.end(); ++iter) {
	if (i == *iter-1){
	  _useLayersBarrelVec[i]=true; break;
	}
      }
    }
  }

  if(_layersToKeepEndcapVec.size() > 0) {
    //layers start at 0
    for(int i = 0; i < layersEndcap; ++i) {
      _useLayersEndcapVec.push_back(false);
      for(IntVec::iterator iter = _layersToKeepEndcapVec.begin(); iter < _layersToKeepEndcapVec.end(); ++iter) {
	if (i == *iter-1){
	  _useLayersEndcapVec[i]=true; break;
	}
      }
    }
  }
// --- Get the muon detector IDs
  dd4hep::Detector &theDet = dd4hep::Detector::getInstance();
  _muonDetBarrel = theDet.constant<unsigned int>("DetID_Yoke_Barrel");
  _muonDetEndcap = theDet.constant<unsigned int>("DetID_Yoke_Endcap");

  _rng = gsl_rng_alloc(gsl_rng_ranlxs2);
  marlin::Global::EVENTSEEDER->registerProcessor(this);
}


void DDSimpleMuonDigi::processRunHeader( LCRunHeader* /*run*/) {
  _nRun++ ;
  _nEvt = 0;
} 

void DDSimpleMuonDigi::processEvent( LCEvent * evt ) { 
    
  // --- Seed the random generator engine
  gsl_rng_set( _rng, marlin::Global::EVENTSEEDER->getSeed(this) ) ;

  streamlog_out( DEBUG ) << " process event : " << evt->getEventNumber() 
			 << " - run  " << evt->getRunNumber() << std::endl ;


  LCCollectionVec *muoncol = new LCCollectionVec(LCIO::CALORIMETERHIT);
  // Relation collection CalorimeterHit, SimCalorimeterHit
  LCCollection* chschcol = 0;
  UTIL::LCRelationNavigator calohitNav = UTIL::LCRelationNavigator( LCIO::CALORIMETERHIT, LCIO::SIMCALORIMETERHIT );

  LCFlagImpl flag;

  flag.setBit(LCIO::CHBIT_LONG);
  flag.setBit(LCIO::CHBIT_ID1);
  flag.setBit(LCIO::RCHBIT_TIME);

  muoncol->setFlag(flag.getFlag());

  // 
  // * Reading Collections of MUON Simulated Hits * 
  // 
  string initString;
  for (unsigned int i(0); i < _muonCollections.size(); ++i) {

    std::string colName =  _muonCollections[i] ;
    
    //fg: need to establish the yoke subdetetcor part here 
    //    use collection name as cellID does not seem to have that information
    CHT::Layout caloLayout = layoutFromString( colName ) ; 

    try{

      LCCollection * col = evt->getCollection( _muonCollections[i].c_str() ) ;
      initString = col->getParameters().getStringVal(LCIO::CellIDEncoding);
      int numElements = col->getNumberOfElements();
      CellIDDecoder<SimCalorimeterHit> idDecoder( col );

      for (int j(0); j < numElements; ++j) {

	      SimCalorimeterHit * hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt( j ) ) ;
	      float energy = hit->getEnergy();
        int cellid = hit->getCellID0();
        int cellid1 = hit->getCellID1();
        unsigned int system = idDecoder(hit)["system"];

        //Get The LayerNumber 
        unsigned int layer = abs( idDecoder(hit)[ _cellIDLayerString ] ) ;
        //Check if we want to use this layer, else go to the next hit
        if( !useLayer(caloLayout, layer) ) continue;

        //calibrate energy deposition
        float calibr_coeff(1.);
        calibr_coeff = _calibrCoeffMuon;
        float hitEnergy = calibr_coeff*energy;
        if(hitEnergy>_maxHitEnergyMuon)hitEnergy=_maxHitEnergyMuon;

        //get hit time, smear and check against time window
        float hitTime = computeHitTime(hit);
        hitTime += gsl_ran_gaussian(_rng, _timeResolution);

        bool hitInTimeWindow = timeHitCut(hit->getPosition(), hitTime);

        //fill if hit passes energy and time filter
        if (hitEnergy > _thresholdMuon && hitInTimeWindow) {
          CalorimeterHitImpl * calhit = new CalorimeterHitImpl();
          calhit->setCellID0(cellid);
          calhit->setCellID1(cellid1);
          calhit->setEnergy(hitEnergy);
          calhit->setPosition(hit->getPosition());
          calhit->setType( CHT( CHT::muon, CHT::yoke, caloLayout ,  idDecoder(hit)[ _cellIDLayerString ] ) );
          calhit->setTime(hitTime);
          calhit->setRawHit(hit);
          muoncol->addElement(calhit);
          calohitNav.addRelation(calhit, hit, 1.0);
        }

      }
    }
    catch(DataNotAvailableException &e){ 
    }
  }
  muoncol->parameters().setValue(LCIO::CellIDEncoding,initString);
  evt->addCollection(muoncol,_outputMuonCollection.c_str());
  // Create and add relation collection for ECAL/HCAL to event
  chschcol = calohitNav.createLCCollection();
  evt->addCollection(chschcol,_outputRelCollection.c_str());


  _nEvt++;

}


void DDSimpleMuonDigi::check( LCEvent * /*evt*/ ) { }
  
void DDSimpleMuonDigi::end(){ } 

bool DDSimpleMuonDigi::useLayer(CHT::Layout caloLayout,  unsigned int layer) {
  switch (caloLayout){
  case CHT::barrel:
    if(layer > _useLayersBarrelVec.size() || _useLayersBarrelVec.size() == 0) return true;
    return _useLayersBarrelVec[layer]; //break not needed, because of return
  case CHT::endcap:
    if(layer > _useLayersEndcapVec.size() || _useLayersEndcapVec.size() == 0) return true;
    return _useLayersEndcapVec[layer]; //break not needed, because of return
  //For all other cases, always keep the hit
  default:
    return true;
  }
}//useLayer


float DDSimpleMuonDigi::computeHitTime( const EVENT::SimCalorimeterHit *h ) const {
  if( nullptr == h ) {
    return 0.f ;
  }
  // Sort sim hit MC contribution by time.
  // Accumulate the energy from earliest time till the energy
  // threshold is reached. The hit time is then estimated at this position in the array
  const float* pos = h->getPosition();
  float deltaT = std::sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]) / _cmm;
  
  using entry_type = std::pair<float, float> ;
  std::vector<entry_type> timeToEnergyMapping {} ;

  const unsigned int nContribs = h->getNMCContributions() ;
  timeToEnergyMapping.reserve( nContribs ) ;

  for( unsigned int c=0 ; c<nContribs ; ++c ) {
    timeToEnergyMapping.push_back( { h->getTimeCont(c), h->getEnergyCont(c) } ) ;
  }
  std::sort(timeToEnergyMapping.begin(), timeToEnergyMapping.end(), [this](entry_type &lhs, entry_type &rhs){
      return lhs.first < rhs.first ;
  }) ;
  float energySum = 0.f ;
  for(auto &entry : timeToEnergyMapping ) {
    energySum += entry.second * _calibrCoeffMuon;
    if( energySum > _timeThresholdMuon ) {
      return entry.first - deltaT; //correct for ToF
    }
  }
  
  return 999.f ; //threshold not overcome
}

// void DDSimpleMuonDigi::smearPosition(const float* pos, float* corr_pos){

//   //cylindrical coordinate sys
//   float x = pos[0];
//   float y = pos[1];
//   float z = pos[2];
//   float r = std::sqrt(x*x+y*y);
//   float phi = std::acos(x/r); 

//   float phiResol = _xyResolution / r;

//   //smear the zPhi, assuming the det is approx a cylinder
//   float phiSmear = phi + gsl_ran_gaussian(_rng, phiResol);
//   float zSmear = z + gsl_ran_gaussian(_rng, _xyResolution);

//   float xSmear = r * std::cos(phiSmear);
//   float ySmear = (y > 0. ? 1. : -1.) * r * std::sin(phiSmear);

//   corr_pos[0] = xSmear;
//   corr_pos[1] = ySmear;
//   corr_pos[2] = zSmear;

//   return;
// }

bool DDSimpleMuonDigi::timeHitCut(const float* pos, float t){

  float x = pos[0], y = pos[1], z = pos[2];
  float theta = std::acos( z / std::sqrt( x*x + y*y + z*z ) ) * 180. / 3.14159;
  theta = theta < 90. ? theta : 180. - theta;

  int idx = 0;
  for(int i = 0; i < _angleRegions.size(); i++){
    if(theta < _angleRegions[i]) break;
    idx++;
  }

  if(t > _timeMins[idx] && t < _timeMaxs[idx]) return true;
  else return false;
}