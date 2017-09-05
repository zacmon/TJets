/*----------------------------------------------------------------------

TITLE : MakeNTupleFromPythia.cc

DESCRIPTION : This runs the Pythia8 event generator, parsing the particle list
to determine which particles are to be used as final state particles [status = 81-99]
and writes these out as a vector of doubles for the components of the 4-vectors.  It further
writes out 4-vectors for each of the truth particles coming from the process itself (partons, W/Z, top)
depending on the process being simulated that can be used for truth identification of jets after clustering

COMPILE :
$ source compile.sh

RUN :
$ ./MakeNTupleFromPythia <type> <output>

type   : 0 = dijet , 1 = ttbar , 2 = W'->WZ (this determines which config file is used)
output : Anything you want - but being logical.  This will be the input for the NTupler code

//----------------------------------------------------------------------*/

#include "MakeNTupleFromPythia.h"
using namespace Pythia8;

int main(int argc, char* argv[]){



  //exit if you dont pass a run card
  if(argc<4){
    cout<<"You need to specify more arguments"<<endl;
    cout<<"Arg1 = process type"<<endl;
    cout<<"Arg2 = output file path and name"<<endl;
    cout<<"Arg3 = number of events"<<endl;
    cout<<"Arg4 = debug flag (optional)"<<endl;
    return 1;
  }

  //get processtype
  std::string ProcessType = argv[1];

  //get config file
  string ConfigFile   = "config_DEFAULT.cmnd";

  if(ProcessType=="dijet")
    ConfigFile="config_dijet.cmnd";
  else if(ProcessType=="ww")
    ConfigFile="config_ww.cmnd";
  else if(ProcessType=="tt")
    ConfigFile="config_ttbar.cmnd";
  else if(ProcessType=="hh")
    ConfigFile="config_higgs.cmnd";
  else{
    cout<<"Bad process type!"<<endl;
    return 0;
  }

  //get output file name
  string OutputFile = argv[2];

  //number of events
  int nEvents = atoi(argv[3]);

  //debug flag
  bool debug=false;
  if(argc>=5){
    string argdebug = argv[4];
    if(argdebug=="debug")
      debug=true;
  }
  //print out the input arguments
  cout<<"InputArguments:  ProcessType="<<ProcessType<<"   ConfigFile="<<ConfigFile<<"  OutputFile="<<OutputFile<<"  Debug="<<debug<<endl;


  ///////////////////////////////////
  //initialize the output ttree and branches
  ///////////////////////////////////


  TTree *tree = new TTree("tree","tree");
  tree->Branch("fspart_id", &fspart_id);
  tree->Branch("fspart_pt", &fspart_pt);
  tree->Branch("fspart_eta",&fspart_eta);
  tree->Branch("fspart_phi",&fspart_phi);
  tree->Branch("fspart_m",  &fspart_m);

  tree->Branch("truth_q1_pt",  &truth_q1_pt);
  tree->Branch("truth_q1_eta", &truth_q1_eta);
  tree->Branch("truth_q1_phi", &truth_q1_phi);
  tree->Branch("truth_q1_m",   &truth_q1_m);

  tree->Branch("truth_q2_pt",  &truth_q2_pt);
  tree->Branch("truth_q2_eta", &truth_q2_eta);
  tree->Branch("truth_q2_phi", &truth_q2_phi);
  tree->Branch("truth_q2_m",   &truth_q2_m);

  tree->Branch("truth_t1_pt",  &truth_t1_pt);
  tree->Branch("truth_t1_eta", &truth_t1_eta);
  tree->Branch("truth_t1_phi", &truth_t1_phi);
  tree->Branch("truth_t1_m",   &truth_t1_m);

  tree->Branch("truth_t2_pt",  &truth_t2_pt);
  tree->Branch("truth_t2_eta", &truth_t2_eta);
  tree->Branch("truth_t2_phi", &truth_t2_phi);
  tree->Branch("truth_t2_m",   &truth_t2_m);

  tree->Branch("truth_W1_pt",  &truth_W1_pt);
  tree->Branch("truth_W1_eta", &truth_W1_eta);
  tree->Branch("truth_W1_phi", &truth_W1_phi);
  tree->Branch("truth_W1_m",   &truth_W1_m);

  tree->Branch("truth_W2_pt",  &truth_W2_pt);
  tree->Branch("truth_W2_eta", &truth_W2_eta);
  tree->Branch("truth_W2_phi", &truth_W2_phi);
  tree->Branch("truth_W2_m",   &truth_W2_m);

  tree->Branch("truth_H_pt",  &truth_H_pt);
  tree->Branch("truth_H_eta", &truth_H_eta);
  tree->Branch("truth_H_phi", &truth_H_phi);
  tree->Branch("truth_H_m",   &truth_H_m);

  //initialize pythia
  Pythia pythia;

  //read in config card
  pythia.readFile(ConfigFile);

  //get the number of events from the config file
  //int nEvents = pythia.mode("Main:numberOfEvents");

  //collide protons with setup defined above
  pythia.init();

  int nAcceptedEvents=0;

  ///////////////////////////////////
  //main event loop
  ///////////////////////////////////
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {

    //print out to show progress
    if(iEvent%100==0)
      cout<<"Pythia: ProcessType="<<ProcessType<<"  entry="<<iEvent<<endl;

    //generates a new event
    pythia.next();

    //accept event flag
    bool acceptevent=false;

    //clear all branches in preparation for next event
    fspart_id.clear();
    fspart_pt.clear();
    fspart_eta.clear();
    fspart_phi.clear();
    fspart_m.clear();

    truth_q1_pt  = 100.0;
    truth_q1_eta = 10.0;
    truth_q1_phi = 0.0;
    truth_q1_m   = -10.0;

    truth_q2_pt  = 100.0;
    truth_q2_eta = 10.0;
    truth_q2_phi = 0.0;
    truth_q2_m   = -10.0;

    truth_t1_pt  = 100.0;
    truth_t1_eta = 10.0;
    truth_t1_phi = 0.0;
    truth_t1_m   = 1.0;

    truth_t2_pt  = 100.0;
    truth_t2_eta = 10.0;
    truth_t2_phi = 0.0;
    truth_t2_m   = 1.0;

    truth_W1_pt  = 100.0;
    truth_W1_eta = 10.0;
    truth_W1_phi = 0.0;
    truth_W1_m   = 1.0;

    truth_W2_pt  = 100.0;
    truth_W2_eta = 10.0;
    truth_W2_phi = 0.0;
    truth_W2_m   = 1.0;

    truth_H_pt  = 100.0;
    truth_H_eta = 10.0;
    truth_H_phi = 0.0;
    truth_H_m   = 1.0;

    bool flag_q1_set=false;

    //loops through the particles in the event just generated
    for (int iPart = 0; iPart < pythia.event.size(); ++iPart) {

      //event filter for
      if(debug) cout<<"Process specific filters: "<<ProcessType<<endl;
      if(ProcessType=="dijet"){
        if(debug) cout<<"Dijet Filter"<<endl;
        //fill truth quark branches
        if(debug) cout<<"FilterStatus : "<<pythia.event[iPart].status()<<"  "<<flag_q1_set<<endl;
        if(pythia.event[iPart].status()==-23 && !flag_q1_set){
          if(debug) cout<<"GotTruth - q1"<<endl;
          truth_q1_pt=pythia.event[iPart].pT();
          truth_q1_eta=pythia.event[iPart].eta();
          truth_q1_phi=pythia.event[iPart].phi();
          truth_q1_m=pythia.event[iPart].m();
          if(debug) cout<<"q1  "<<truth_q1_pt<<"  "<<truth_q1_eta<<"  "<<truth_q1_phi<<"  "<<truth_q1_m<<endl;
          if(truth_q1_pt>190){
            acceptevent=true;
          }

          flag_q1_set=true;

        }
        else if(pythia.event[iPart].status()==-23){
          if(debug) cout<<"GotTruth - q2"<<endl;
          truth_q2_pt=pythia.event[iPart].pT();
          truth_q2_eta=pythia.event[iPart].eta();
          truth_q2_phi=pythia.event[iPart].phi();
          truth_q2_m=pythia.event[iPart].m();
          if(debug) cout<<"q2  "<<truth_q2_pt<<"  "<<truth_q2_eta<<"  "<<truth_q2_phi<<"  "<<truth_q2_m<<endl;
          if(truth_q2_pt>190){
            acceptevent=true;
          }
        }
      }
      else if(ProcessType=="ww"){
        if(debug) cout<<"WW Filter"<<endl;
        //only accept if Wprime decay is via G*->WW - identify by asking for Z0 in the intermediate state
        if(pythia.event[iPart].id()==24 && pythia.event[iPart].status()==-22){
          acceptevent=true;
        }
        //fill truth boson branches
        if(pythia.event[iPart].id()==24 && pythia.event[iPart].status()==-22){
          if(debug) cout<<"GotTruth - W1"<<endl;
          truth_W1_pt  = pythia.event[iPart].pT();
          truth_W1_eta = pythia.event[iPart].eta();
          truth_W1_phi = pythia.event[iPart].phi();
          truth_W1_m   = pythia.event[iPart].m();
        }
        else if( pythia.event[iPart].id()==-24 && pythia.event[iPart].status()==-22){
          if(debug) cout<<"GotTruth - W2"<<endl;
          truth_W2_pt  = pythia.event[iPart].pT();
          truth_W2_eta = pythia.event[iPart].eta();
          truth_W2_phi = pythia.event[iPart].phi();
          truth_W2_m   = pythia.event[iPart].m();
        }
      }
      else if(ProcessType=="tt"){
        if(debug) cout<<"Top quark Filter"<<endl;
        //fill truth topquark branches
        if(pythia.event[iPart].id()==6 && pythia.event[iPart].status()==-62){
          if(debug) cout<<"GotTruth - t1"<<endl;
          truth_t1_pt=pythia.event[iPart].pT();
          truth_t1_eta=pythia.event[iPart].eta();
          truth_t1_phi=pythia.event[iPart].phi();
          truth_t1_m=pythia.event[iPart].m();
          if(truth_t1_pt>180){
            acceptevent=true;
          }
        }
        else if(pythia.event[iPart].id()==-6 && pythia.event[iPart].status()==-62){
          if(debug) cout<<"GotTruth - t2"<<endl;
          truth_t2_pt=pythia.event[iPart].pT();
          truth_t2_eta=pythia.event[iPart].eta();
          truth_t2_phi=pythia.event[iPart].phi();
          truth_t2_m=pythia.event[iPart].m();
          if(truth_t2_pt>180){
            acceptevent=true;
          }
        }
      }
      else if(ProcessType=="hh"){
        if(debug) cout<<"Higgs Filter"<<endl;
        //fill truth higgs branches
        if(pythia.event[iPart].id()==25 && pythia.event[iPart].status()==-62){
          if(debug) cout<<"GotTruth - H"<<endl;
          truth_H_pt=pythia.event[iPart].pT();
          truth_H_eta=pythia.event[iPart].eta();
          truth_H_phi=pythia.event[iPart].phi();
          truth_H_m=pythia.event[iPart].m();
          if(truth_H_pt>180){
            acceptevent=true;
          }
        }
      }
      else{
        if(debug) cout<<"Filling event"<<endl;
        acceptevent=true;
      }
      if(debug) cout<<acceptevent<<endl;


      //only save particles to output ttree if they are final state particles
      //this is the method recommended by Yang-Ting
      if(pythia.event[iPart].isFinal()){
        if(debug) {
          if(iEvent%10==0){
            cout<<pythia.event[iPart].status()<<"  "
                <<pythia.event[iPart].id()<<"  "
                <<pythia.event[iPart].pT()<<"  "
                <<pythia.event[iPart].eta()<<"  "
                <<pythia.event[iPart].phi()<<"  "
                <<pythia.event[iPart].m()<<endl;
          }
        }

        fspart_id .push_back(pythia.event[iPart].id());
        fspart_pt .push_back(pythia.event[iPart].pT());
        fspart_eta.push_back(pythia.event[iPart].eta());
        fspart_phi.push_back(pythia.event[iPart].phi());
        fspart_m  .push_back(pythia.event[iPart].m());
      }
    }

    //fill event into tree if it has passed the event filter
    if(acceptevent){
      nAcceptedEvents+=1;
      tree->Fill();
    }


  }

  ///////////////////////////////////
  //Summary of truth filters
  ///////////////////////////////////
  cout<<"ProcessedEvents:"<<endl
      <<"total      = "<<nEvents<<endl
      <<"accepted   = "<<nAcceptedEvents<<endl
      <<"efficiency = "<<float(nAcceptedEvents)/float(nEvents)<<endl;

  TFile *fout = new TFile(OutputFile.c_str(),"RECREATE");
  tree->Write("tree");
  fout->Close();

  return 0;
}
