#include "gobbi.h"

//Written by Nicolas Dronchi
//Modified by Johnathan Phillips 2024_08_29
//Class written to handle all specifics of the gobbi array
//such as communicating with HINP, calibrations, Checking for charge sharing in neighbor
//calculating geometry
//                 ,,__
//        ..  ..   / o._)                   .---.
//       /--'/--\  \-'||        .----.    .'     '.
//      /        \_/ / |      .'      '..'         '-.
//    .'\  \__\  __.'.'     .'          e-._
//      )\ |  )\ |      _.'
//     // \\ // \\
//    ||_  \\|_  \\_             (Wait Gobbi the person, not Gobi desert)
//    '--' '--'' '--'


//TODO
//     Only want matching for F, B, and CsI
//     All Eloss needs to be for new setup, have aluminum absorbers now

//TODO write an unpacker inside of Gobbi that calls HINP for Si and ADC and TDC for CsI
//The unpacker should accept a source ID which will tell it if it goes to HINP, ADC, or TDC

//VME crate order
//	Timestamper (SIS)
//	XLM (HINP boards)
//	CAEN ADC (CsI)
//	CAEN QDC (CsI)
//	CAEN TDC (CsI)
//	CAESAR (Called separately in det.cpp)


using namespace std;

gobbi::gobbi(TRandom* ran, histo_sort * Histo1)
{
  Histo = Histo1;
                        //Ntele, Nstrip, filename, polynomial order, weave strips?

                       //Ntele0, Nstrip0,      name,                 order0, weave,  bback
  FrontEcal = new calibrate(4, Histo->Nstrip, "cal2/GobbiFrontEcal.dat", 1, false);
  BackEcal = new calibrate(4, Histo->Nstrip, "cal2/GobbiBackEcal.dat", 1, false);
  //CsIEcal = new calibrate(4, 4, "cal/calMacros/CsI_cal_out_0deg.dat", 1, false); //0 degree energies
  CsIEcal = new calibrate(4, 4, "cal/calMacros/CsI_cal_out.dat", 1, false); //Center of CsI angle
  //CsIEcal = new calibrate(4, 4, "cal/CsI_roughcal.dat", 1, false); //Rough calibration using the pedestal and the highest energy proton point
  FrontTimecal = new calibrate(4, Histo->Nstrip, "cal/GobbiFrontTimecal.dat",1, false);
  BackTimecal = new calibrate(4, Histo->Nstrip, "cal/GobbiBackTimecal.dat",1, false);  
  CsITimecal = new calibrate(4, 4, "cal/GobbiCsITimecal.dat", 1, false);

  FrontLowEcal = new calibrate(4, Histo->Nstrip, "cal2/GobbiFrontLowLowEcal.dat", 1, false);
  BackLowEcal = new calibrate(4, Histo->Nstrip, "cal2/GobbiBackLowEcal.dat", 1, false);

  for (int id=0;id<4;id++)
  {
    Telescope[id] = new telescope(false); //false for not S800
    Telescope[id]->init(id); //tells Telescope what position it is in
    Telescope[id]->load(0,0,15,15,   
                        16,16,31,31,
                        0,15,15,0,  //These values give slightly fuzzy quadrants, intentional
                        16,31,31,16);  //load strip limts for each CsI crystal
  }

  //TODO don't think I need to do init() for this, but check that
  S800 = new telescope(false); //Create 5th "telescope" to hold S800 solutions, true for S800
  S800->init(0);

  
  for (int i=0;i<4;i++)
  {
    Gfvb_cnt[i]=0; //tracks f vs b for plotting
  }
  
  SiADC = new HINP();
  CsIADC = new caen();
  CsITDC = new mtdc();
  CsIQDC = new caen();

  //time gates for s800 coincidence
  int one,two,three;
  ifstream fileTime("cal/csiTime.gate");
  for (int i=0;i<20;i++)
  {
    fileTime >> one >> two >> three;
    csiTimeMin[one]= two;
    csiTimeMax[one]= three;
  }

}

gobbi::~gobbi()
{
  delete FrontEcal;
  delete BackEcal;
  delete CsIEcal;
  delete FrontTimecal;
  delete BackTimecal;
  delete CsITimecal;
  delete FrontLowEcal;
  delete BackLowEcal;
  //delete [4] Telescope; //not needed as it is in automatic memory, didn't call with new
  delete S800; //not needed as it is in automatic memory, didn't call with new
  delete CsITDC;
  delete CsIADC;
  delete CsIQDC;
  delete SiADC; //added in these deletes, not sure they should be here?
}

void gobbi::SetTarget(double Targetdist, float TargetThickness)
{
  //TODO Target position will need to change with S800 setting
  for (int id=0;id<4;id++){ Telescope[id]->SetTarget(Targetdist, TargetThickness); }

  //S800 is how far from the target? Will be const
  //Guess 70 cm?, TODO must get actual distance
  S800->SetTarget(70., TargetThickness);
}
void gobbi::reset()
{
  //reset the Telescope class
  for (int i=0;i<4;i++){ Telescope[i]->reset();}
  S800->reset();
  //make sure to reset the CsI here as well, whatever they end up looking like
}

//TODO connector unpacking class
//The merged file will have one source ID for the VME crate, need to unpack HINP, CsI ADC, and CsI TDC
//unpack() calls sub-unpackers for each of these
//Not sure what the file structure looks like between these different VME modules.
// ****** Right now we have to comment out HINP and TDC, don't want to do that ******

//TODO do I actually need tdcpoint? Remove if possible
bool gobbi::unpack(unsigned short *&point,unsigned short *&tdcpoint,int runno)
{

  //TODO fabricating qdc data, hardcode zeros at start for qdc flag, seems to work
  for (int i=0;i<32;i++) DataE[i].qdcmatch = 0;

  NE = 0;
  NT = 0;
  NQ = 0; //hardcode these as 0, makes it much safer

  bool stat = true;
  stat = unpackSi_HINP4(point);
  if (!stat)
  {
    cout << "Bad Si-HINP" << endl;
    return stat;
  }

  //TODO is a shift needed here?
  //point += 3;
  //cout << "NEW RUN" << endl;
  stat = unpackCsI_ADC(point);
 
  if (!stat)
  {
    //cout << "Bad CsI ADC" << endl;
    return stat;
  }
	stat = unpackCsI_QDC(point); //QDC basically uses the same unpacker as ADC

	if (!stat)
	{
		cout << "Bad CsI QDC" << endl;
		return stat;
	}

  stat = unpackCsI_TDC(point);

  if (!stat)
  {
    cout << "Bad CsI TDC" << endl;
    return stat;
  }

	//Calculate CsI PSD parameter first, then match with time
	CsIQDCMatch();
  //Match up CsI energies and times now
  MatchCsIEnergyTime();

  return stat;
}

//TODO unpacks the silicon HINP boards, 8 used
bool gobbi::unpackSi_HINP4(unsigned short*& point)
{
  bool stat = true;
  stat = SiADC->unpackSi_HINP4(point); //unpack all the HINP chips
  if (!stat)
  {
    cout << "Bad Si-HINP" << endl;
    //cntBadHinp++;
    return false;
  }

  //TODO After unpacking, add front and back events
  for (int i=0; i<SiADC->NstripsRead; i++)
  {
    float Energy = 0;
    float time = 0; //can be calibrated or shifted later

    //position in board tells you what you have and what quadrant you are in
    //Add time gates once you get an idea of the timing distributions
    //All Gobbi SiE detectors
    if (SiADC->board[i] == 1 || SiADC->board[i] == 3 || SiADC->board[i] == 5 || SiADC->board[i] == 7)
    {
      int quad = (SiADC->board[i] - 1)/2; 
      addFrontEvent(quad, SiADC->chan[i], SiADC->high[i], SiADC->low[i], SiADC->time[i]);
    }
    if (SiADC->board[i] == 2 || SiADC->board[i] == 4 || SiADC->board[i] == 6 || SiADC->board[i] == 8)
    {
      int quad = (SiADC->board[i]/2)-1;
      addBackEvent(quad, SiADC->chan[i], SiADC->high[i], SiADC->low[i], SiADC->time[i]);
    }
  }

  return stat;

}

//TODO unpacks the CsI ADC channels, 16 used
bool gobbi::unpackCsI_ADC(unsigned short*& point)
{
  NE = 0; //index for the current energy position

  bool stat = true;
  //look to see if there is any CAEN ADC data at all
  unsigned short header1 = *point;
  unsigned short header2 = *(point+1);  
	//cout << hex << header1 << " " << header2 << dec << endl;
  if (header1 == 0xffff)
  {
    //cout << " no data " << endl;
    CsIADC->number = 0;   //  no data at all
    //cntNoAdc++;
    *point+=2;
    return false;
  }
  else  // read data
  {
    point = CsIADC->read(point); //unpack data from caen ADC
    if (!stat)
    {
      cout << "Bad CsIE " << CsIADC->number << " " << counter << endl;
      //cntBadAdc++;
      return stat;
    }
  }

  //TODO add CsI events
  if (Verb) {cout << "CsIADC->number " << CsIADC->number << "TDC->Ndata " << CsITDC->number << endl;}


  //cout << " ADC Num " << CsIADC->number << endl;
  for (int i=0; i<CsIADC->number; i++)
  {
    if (CsIADC->underflow[i]) cout << "underflow!" << endl;//continue;
    if (CsIADC->overflow[i]) cout << "overflow!" << endl;// continue;
    
    if(CsIADC->channel[i] < 16) //CsI used chan 16-31 in the ADC
    {
      //want to shift it to 0-15 for it's use in code
      int idtemp = CsIADC->channel[i];
      DataE[NE].id = idtemp;
      DataE[NE].ienergy = CsIADC->data[i];
      NE++;
    }

  }

  unsigned short f1 = *point++;
  unsigned short f2 = *point++;  
	//cout << "ADC end " << hex << header1 << " " << header2 << " ";
  if (f1 != 0xffff && f2 != 0xffff)
  {
		cout << "ADC not read correctly" << endl;
		return false;
  }
  //cout << " end adc " << *point << " " << *(point+1) << endl;

  return stat;
}

//TODO unpacks the CsI TDC channels, 16 used
bool gobbi::unpackCsI_TDC(unsigned short*& point)
{
  NT = 0; //index for the current time position

  bool stat = true;
  CsITDC->number = 0;  //temporarily
  //TODO it is unclear, I think the TDC point is at the end after the CsI ADC????

  //look for FFFF's
  unsigned short header1 = *point;
  unsigned short header2 = *(point + 1);
  //cout << "TDC header " << header1 << " " << header2 << endl;
  if (header1  == 0xffff || header2 == 0xffff)
  {
    //cout << "missing FFFF " << endl;
    //return false;
    //cntNoTDCblock++;
  }
  else if (Verb) cout << "good we found the FFFF blocks" << endl;

  
  unsigned short *checkpoint1 = point;
  unsigned short f1 = *checkpoint1++;
  unsigned short f2 = *checkpoint1++;
  if (f1 != 0xffff && f2 != 0xffff)
  {
    point = CsITDC->read(point); //unpack the mesytec TDC

    unsigned short *checkpoint2 = point;
    f1 = *checkpoint2++;
    f2 = *checkpoint2++;
    if (f1 != 0xffff && f2 != 0xffff)
    {
      //cntNoTDCend++;
      //itHappened++;
      //if (itHappened == 3) abort();
      //cout << "did not read the CsI TDC correctly" << endl;
    }
  }
  else
  {
    //cout << "no CsI TDC data" << endl;
    //point += 2; +=2 might not work?
    point++;
    point++;
    //cntNoTDCdata++;
    return false;
  }

  //TODO add TDC events

  for (int i=0; i<CsITDC->number; i++) //TODO check this loop
  {
    //if(Mtdc[itdc].channel[j] == 34) continue; //Extended timestamp
    if(CsITDC->channel[i] == 32 || CsITDC->channel[i] == 33) { //Trigger Channel
	  CsITDC->trigger = true;
	  continue; }
    if (Verb) {cout << "TDC in det.cpp " << CsITDC->channel[i] << " " << CsITDC->data[i] << endl;}
      //if (CsITDC->channel[i] < 16 && CsITDC->channel[i] > 3) continue;
      if (CsITDC->channel[i] < 16 && CsITDC->channel[i] != 1 && CsITDC->channel[i] != 2 && CsITDC->channel[i] != 3) continue;
      int id = CsITDC->channel[i];
      int itime =CsITDC->data[i];
      //if (id == 16) cout << "chan 0 time " << endl;
      if (CsITDC->channel[i] > 15) id = id - 16; //CsI
      if (CsITDC->channel[i] == 1) id = 16; //S800 DBS
      if (CsITDC->channel[i] == 2) id = 17; //S800 Object
      if (CsITDC->channel[i] == 3) id = 18; //S800 RF
      DataT[NT].id  = id;
      DataT[NT].itime = itime;
      DataT[NT].time = itime * CsITDC->res / 1000.0; // Convert to ns
      DataT[NT].ichan = CsITDC->channel[i];
      //want to shift it to 0-15 for it's use in code
      NT++;
  }
  if (Verb) {cout << "NE " << NE << " NT " << NT << endl;}

  f1 = *point++;
  f2 = *point++;  
	//cout << header1 << " " << header2 << " ";
  if (f1 != 0xffff && f2 != 0xffff)
  {
		cout << "ADC not read correctly" << endl;
		return false;
  }

  return stat;
}

//TODO unpacks the CsI QDC channels, 16 used. Same unpacker as ADC
bool gobbi::unpackCsI_QDC(unsigned short*& point)
{
  bool stat = true;

  //Reset QDC
  for (int i=0;i<32;i++) DataQ[i].ienergy = 0;

  NQ = 0; //index for the current energy position
  //look to see if there is any CAEN QDC data at all
  unsigned short header1 = *point;
  unsigned short header2 = *(point+1);  

	//cout << "QDC Startheader " << hex << " " << header1 << dec << endl;
  if (header1 == 0xffff)
  {
    //cout << " QDC no data " << header1 << " " << header2 << endl;
    CsIQDC->number = 0;   //  no data at all
    //cout << "p1 " << *point << "p2 " << *(point+1) << "p3 " << *(point+2) << endl;
    *point++;
    *point++;
    //cout << "p2 " << *point << endl;
    //cntNoAdc++;
    return true;
  }
  else  // read data
  {
    point = CsIQDC->read(point); //unpack data from caen ADC
    if (!stat)
    {
      cout << "Bad CsI QDC " << CsIQDC->number << " " << counter << endl;
      //cntBadAdc++;
      return false;
    }
  }

  //TODO add CsI events
  if (Verb) {cout << "CsIQDC->number " << CsIQDC->number << endl;}

  //cout << endl << CsIQDC->number << endl;
  for (int i=0; i<CsIQDC->number; i++)
  {
    if (CsIQDC->underflow[i]) continue;
    if (CsIQDC->overflow[i]) continue;
    
    if(CsIQDC->channel[i] < 17) //CsI used chan 0-15 in the QDC chan 17 for S800 testing
    {
      //want to shift it to 0-15 for it's use in code
      int idtemp = CsIQDC->channel[i];
      DataQ[NQ].id = idtemp;
      DataQ[NQ].ienergy = CsIQDC->data[i];
      NQ++;

    }

  }


  unsigned short f1 = *point++;
  unsigned short f2 = *point++;  
	//cout << "QDC end " <<  header1 << " " << header2 << " ";
  if (f1 != 0xffff && f2 != 0xffff)
  {
		cout << "QDC not read correctly" << endl;
		return false;
  }

  return stat;
}

void gobbi::addFrontEvent(int quad, unsigned short chan, unsigned short high, 
                                    unsigned short low, unsigned short timeR)
{

  //Use calibration to get Energy and fill elist class in Telescope
  float Energy = FrontEcal->getEnergy(quad, chan, high);
  float EnergyLow = FrontLowEcal->getEnergy(quad, chan, low);
  float Vlow = 0;

  float time = FrontTimecal->getTime(quad, chan, timeR);

  //good spot to fill in all of the summary and chan spectrums
  Histo->sumFrontE_R->Fill(quad*Histo->Nstrip + chan, high);
  Histo->sumFrontlowE_R->Fill(quad*Histo->Nstrip + chan, low);
  Histo->sumFrontTime_R->Fill(quad*Histo->Nstrip + chan, timeR);

  Histo->sumFrontE_cal->Fill(quad*Histo->Nstrip + chan, Energy);
  Histo->sumFrontlowE_cal->Fill(quad*Histo->Nstrip + chan, EnergyLow);
  Histo->sumFrontTime_cal->Fill(quad*Histo->Nstrip + chan, time);
  //if (quad == 0 && chan == 0 && high > 1270 && high < 1340) abort();
  Histo->FrontE_R[quad][chan]->Fill(high);
  Histo->FrontElow_R[quad][chan]->Fill(low);

  Histo->FrontTime_R[quad][chan]->Fill(timeR);
  Histo->FrontE_cal[quad][chan]->Fill(Energy);   


  //TODO this is a good spot to throw an if statement and make software thresholds
  if (Energy > 0.5)
  {
    Telescope[quad]->Front.Add(chan, Energy, EnergyLow, low, high, time, -1,0); //No psd, set to -1
    Telescope[quad]->tempFront.Add(chan, Energy, EnergyLow, low, high, time, -1,0); //TODO what's this for? Added later in Ne16?
    Telescope[quad]->multFront++;
  }
}

void gobbi::addBackEvent(int quad, unsigned short chan, unsigned short high, 
                                   unsigned short low, unsigned short timeR)
{

  //Use calibration to get Energy and fill elist class in Telescope
  float Energy = BackEcal->getEnergy(quad, chan, high);
  float EnergyLow = BackLowEcal->getEnergy(quad, chan, low);

  float time = BackTimecal->getTime(quad, chan, timeR);

  //good spot to fill in all of the summary and chan spectrums
  Histo->sumBackE_R->Fill(quad*Histo->Nstrip + chan, high);
  Histo->sumBackTime_R->Fill(quad*Histo->Nstrip + chan, timeR);
  Histo->sumBackE_cal->Fill(quad*Histo->Nstrip + chan, Energy);
  Histo->sumBackTime_cal->Fill(quad*Histo->Nstrip + chan, time);

  Histo->BackE_R[quad][chan]->Fill(high);

  Histo->BackElow_R[quad][chan]->Fill(low);
  Histo->BackTime_R[quad][chan]->Fill(timeR);

  Histo->BackE_cal[quad][chan]->Fill(Energy); 
  Histo->sumBacklowE_cal->Fill(quad*Histo->Nstrip + chan, EnergyLow);

  Histo->sumBacklowE_cal->Fill(quad*Histo->Nstrip + chan, EnergyLow);
  Histo->sumBacklowE_R->Fill(quad*Histo->Nstrip + chan, low);

  //TODO this is a good spot to throw an if statement and make software thresholds
  if (Energy > 0.5)
  {
    Telescope[quad]->Back.Add(chan, Energy, EnergyLow, low, high, time, -1,0);
    Telescope[quad]->tempBack.Add(chan, Energy, EnergyLow, low, high, time, -1,0); //TODO what's this for? Added later in Ne16?
    Telescope[quad]->multBack++;
  }
}


void gobbi::MatchCsIEnergyTime() //TODO verify that this code works, must turn stuff off if no TDC
{
  int Nfound = 0;
  int Nnotfound = 0;

  //plot unmatched energy and times
  for (int ie=0;ie<NE;ie++)
  {
    double Ecal = CsIEcal->getEnergy(floor(DataE[ie].id/4), DataE[ie].id - (4*floor(DataE[ie].id/4)), DataE[ie].ienergy);

    Histo->CsI_Energy_R_um[DataE[ie].id]->Fill(DataE[ie].ienergy);
    Histo->CsI_Energy_cal[DataE[ie].id]->Fill(Ecal);
    Histo->sumCsIE_R->Fill(DataE[ie].id, DataE[ie].ienergy);
    Histo->sumCsIE_cal->Fill(DataE[ie].id, Ecal);
    if (DataE[ie].qdcmatch == 1) Histo->CsI_QDC_matched[DataE[ie].id]->Fill(DataE[ie].qdc);
  }
  for (int it=0;it<NT;it++)
  {
    //cout << DataT[it].id << endl;
    Histo->CsI_Time_R_um[DataT[it].id]->Fill(DataT[it].time);
    Histo->sumCsITime_R->Fill(DataT[it].id, DataT[it].itime);
  }
  
  //TODO Temporarily don't match times for alpha calibrations. No TDC data yet. Change once we have TDC data
  for (int ie=0;ie<NE;ie++)
  {
  	int id = DataE[ie].id;
 	if (id <4)
        {
          //ignores events with CsI output lower than pedestal, only keep for quench fit
          //if (id == 0 && DataE[ie].ienergy < CsIEcal->Coeff[0][0].qp) continue;
          int quad = 0;
          int chan = DataE[ie].id;
          //CsI energy calibration here!
          double Ecal = CsIEcal->getEnergy(quad, chan, DataE[ie].ienergy);
          //Ecal = Ecal*1.5; // temporarily increase energy to try and straighten theta_H vs erel, modify the pedestal and high energy edge rough cal
          DataE[ie].energy = Ecal;
          DataE[ie].time = CsITimecal->getTime(quad, chan, DataE[ie].itime);

          Telescope[quad]->CsI.Add(chan, DataE[ie].energy, 0., 0, 
                                   DataE[ie].ienergy, DataE[ie].time, DataE[ie].qdc, DataE[ie].qdcmatch);
          Telescope[quad]->multCsI++;
         //cout << " gobbi:ll " << id << " tele=" << quad << " iCsI=" << chan << endl;
        }
        else if(id < 8 and id >= 4)
        {
          //ignores events with CsI output lower than pedestal
          //if (id == 4 && DataE[ie].ienergy < CsIEcal->Coeff[1][0].qp) continue;
          int quad = 1;
          int chan = DataE[ie].id - 4; //shift id so it is 0,1,2,3
          double Ecal = CsIEcal->getEnergy(quad, chan, DataE[ie].ienergy);
          //Ecal = Ecal*1.5; // temporarily increase energy to try and straighten theta_H vs erel
          DataE[ie].energy = Ecal;
          DataE[ie].time = CsITimecal->getTime(quad, chan, DataE[ie].itime);

          Telescope[quad]->CsI.Add(chan, DataE[ie].energy, 0., 0, 
                                   DataE[ie].ienergy, DataE[ie].time, DataE[ie].qdc, DataE[ie].qdcmatch);
          Telescope[quad]->multCsI++;
           //cout << " gobbi:ll " << id << " tele=" << quad << " iCsI=" << chan << endl;
        }
        else if(id < 12 and id >= 8)
        {
          //ignores events with CsI output lower than pedestal
          //if (id == 8 && DataE[ie].ienergy < CsIEcal->Coeff[2][0].qp) continue;
          int quad = 2;
          int chan = DataE[ie].id - 8; //shift id so it is 0,1,2,3
          double Ecal = CsIEcal->getEnergy(quad, chan, DataE[ie].ienergy);
          //Ecal = Ecal*1.6;
          DataE[ie].energy = Ecal;
 
          DataE[ie].time = CsITimecal->getTime(quad, chan, DataE[ie].itime);

          Telescope[quad]->CsI.Add(chan, DataE[ie].energy, 0., 0, 
                                   DataE[ie].ienergy, DataE[ie].time, DataE[ie].qdc, DataE[ie].qdcmatch);
          Telescope[quad]->multCsI++;
         //cout << " gobbi:ll " << id << " tele=" << quad << " iCsI=" << chan << endl;
        }
        else if(id < 16 and id >= 12)
        {
          //ignores events with CsI output lower than pedestal
          //if (id == 12 && DataE[ie].ienergy < CsIEcal->Coeff[3][0].qp) continue;
          int quad = 3;
          int chan = DataE[ie].id - 12; //shift id so it is 0,1,2,3
          double Ecal = CsIEcal->getEnergy(quad, chan, DataE[ie].ienergy);
          //Ecal = Ecal*1.45; // temporarily increase energy to try and straighten theta_H vs erel
          DataE[ie].energy = Ecal;
          DataE[ie].time = CsITimecal->getTime(quad, chan, DataE[ie].itime);

          Telescope[quad]->CsI.Add(chan, DataE[ie].energy, 0., 0, 
                                   DataE[ie].ienergy, DataE[ie].time, DataE[ie].qdc, DataE[ie].qdcmatch);
          Telescope[quad]->multCsI++;
         //cout << " gobbi:ll " << id << " tele=" << quad << " iCsI=" << chan << endl;
        }
        else
        {
          //we have an issue. what did YOU do?
          cout << "check the id of CsI coming in, CsI.id" << id << " was found" << endl;
        }
  
  }

  /*// match up energies to times
  for (int ie=0;ie<NE;ie++)
  {
    DataE[ie].itime = -1;
    //DataQ[ie].itime = -1;	//This probably doesn't need time.
    bool found = false;
    for (int it=0;it<NT;it++)
    {
      if (found != true && true) //Use this until you get good times. Always matches TDC and CsI
      // TODO We have matched with some time gate imposed - ONLY UNCOMMENT WHEN YOU HAVE GOOD TIME GATES
      //if (DataE[ie].id == DataT[it].id && DataT[it].itime > -599 && DataT[it].itime < 500)
      {
        found = true;
        DataE[ie].itime = DataT[it].itime;
        //DataQ[ie].itime = DataT[it].itime;t

        Nfound++;
        
        //add event to the right telescope elsit of csi 
        //TODO take out energy calibration here, Calibrations need to be for each CsI crystal
        //and is also dependent on PID in CsI. -ND
        int id = DataE[ie].id;

         //test recalibratiopn of CsI energiesroot  ******************* remove eventually
        double fact = 1.2;

        if (id <4)
        {
          //ignores events with CsI output lower than pedestal, only keep for quench fit
          //if (id == 0 && DataE[ie].ienergy < CsIEcal->Coeff[0][0].qp) continue;
          int quad = 0;
          int chan = DataE[ie].id;
          //CsI energy calibration here!
          double Ecal = CsIEcal->getEnergy(quad, chan, DataE[ie].ienergy);

          DataE[ie].energy = Ecal;
          DataE[ie].time = CsITimecal->getTime(quad, chan, DataE[ie].itime);

          Telescope[quad]->CsI.Add(chan, DataE[ie].energy, 0., 0, 
                                   DataE[ie].ienergy, DataE[ie].time, DataE[ie].qdc);
          Telescope[quad]->multCsI++;
         //cout << " gobbi:ll " << id << " tele=" << quad << " iCsI=" << chan << endl;
        }
        else if(id < 8 and id >= 4)
        {
          //ignores events with CsI output lower than pedestal
          //if (id == 4 && DataE[ie].ienergy < CsIEcal->Coeff[1][0].qp) continue;
          int quad = 1;
          int chan = DataE[ie].id - 4; //shift id so it is 0,1,2,3
          double Ecal = CsIEcal->getEnergy(quad, chan, DataE[ie].ienergy);

          DataE[ie].energy = Ecal;
          DataE[ie].time = CsITimecal->getTime(quad, chan, DataE[ie].itime);

          Telescope[quad]->CsI.Add(chan, DataE[ie].energy, 0., 0, 
                                   DataE[ie].ienergy, DataE[ie].time, DataE[ie].qdc);
          Telescope[quad]->multCsI++;
           //cout << " gobbi:ll " << id << " tele=" << quad << " iCsI=" << chan << endl;
        }
        else if(id < 12 and id >= 8)
        {
          //ignores events with CsI output lower than pedestal
          //if (id == 8 && DataE[ie].ienergy < CsIEcal->Coeff[2][0].qp) continue;
          int quad = 2;
          int chan = DataE[ie].id - 8; //shift id so it is 0,1,2,3
          double Ecal = CsIEcal->getEnergy(quad, chan, DataE[ie].ienergy);

          DataE[ie].energy = Ecal;
 
          DataE[ie].time = CsITimecal->getTime(quad, chan, DataE[ie].itime);

          Telescope[quad]->CsI.Add(chan, DataE[ie].energy, 0., 0, 
                                   DataE[ie].ienergy, DataE[ie].time, DataE[ie].qdc);
          Telescope[quad]->multCsI++;
         //cout << " gobbi:ll " << id << " tele=" << quad << " iCsI=" << chan << endl;
        }
        else if(id < 16 and id >= 12)
        {
          //ignores events with CsI output lower than pedestal
          //if (id == 12 && DataE[ie].ienergy < CsIEcal->Coeff[3][0].qp) continue;
          int quad = 3;
          int chan = DataE[ie].id - 12; //shift id so it is 0,1,2,3
          double Ecal = CsIEcal->getEnergy(quad, chan, DataE[ie].ienergy);

          DataE[ie].energy = Ecal;
          DataE[ie].time = CsITimecal->getTime(quad, chan, DataE[ie].itime);

          Telescope[quad]->CsI.Add(chan, DataE[ie].energy, 0., 0, 
                                   DataE[ie].ienergy, DataE[ie].time, DataE[ie].qdc);
          Telescope[quad]->multCsI++;
         //cout << " gobbi:ll " << id << " tele=" << quad << " iCsI=" << chan << endl;
        }
        else
        {
          //we have an issue. what did YOU do?
          cout << "check the id of CsI coming in, CsI.id" << id << " was found" << endl;
        }

        Histo->CsI_Energy_R[id]->Fill(DataE[ie].ienergy);
        Histo->CsI_Energy_cal[id]->Fill(DataE[ie].energy);
        Histo->CsI_Time_R[id]->Fill(DataE[ie].itime);
        Histo->CsI_Time_cal[id]->Fill(DataE[ie].time);

      }
    }

  }*/
}

//Assign the qdc value using the QDC and ADC
void gobbi::CsIQDCMatch()
{
  //cout << endl << "MATCHING EVENTS" << endl;
	for (int i=0;i<NQ;i++)
	{
	  //cout << DataQ[i].id << " ";
	  Histo->CsI_QDC_R[DataQ[i].id]->Fill(DataQ[i].ienergy);
	}
  //cout << endl;
	//channels always read in order
	for (int i=0;i<NE;i++)
	{
    DataE[NE].qdcmatch = 0;
	  //cout << DataE[i].id << " ";
		for (int j=0;j<NQ;j++)
		{
			if (DataE[i].id == DataQ[j].id)
			{
				DataE[i].qdc = DataQ[j].ienergy;
        DataE[i].qdcmatch = 1;
        break;
			}
		}
		//cout << endl;
	}

}


//   _.+._
// (^\/^\/^)
//  \@*@*@/
//  {_____} Just as the queen would prefer to spell it
//Note from 8/29/2024, leave the spelling for old time's sake
void gobbi::SiNeighbours()
{
  //When a charged particle moves through Si, there is cross talk between neighbouring
  //strips. Generally the signal is proportional to the total signal in the Si.
  for (int id=0;id<4;id++) 
  {
    Telescope[id]->Front.Neighbours(id);
    Telescope[id]->Back.Neighbours(id);
  }
}

//TODO need to apply S800 time gates before neighbour adding
//Can't do neighbour adding until I have calibrations
int gobbi::analyze(s800_results S800_results)
{
  //********************************************************

  //TODO Apply S800 time gates here.
  //Do not apply them until needed and with good time gates. Want to accept all solutions for now.

  //********************************************************

  //leave this part out early in the experiment, it causes a lot of trouble!!
  SiNeighbours(); //see if this is working

  //TODO make sure matching works for multiple F, B, CsI. No dE
  //TODO want hists for DEE, F vs B, High vs Low. Think of others to add
  multiECsI = 0;
  int Nmatch = 0;

  for (int id=0;id<4;id++) 
  {
  
    //WARNING: this is only in for testing alpha calibrations. Comment this out
    //when you are taking real data
    //Don't superimpose detectors at different distances with this plot
    if (Telescope[id]->Front.Nstore ==1 && Telescope[id]->Back.Nstore ==1)
    {
      int testhit = Telescope[id]->testingHitE();
      Telescope[id]->position(9);
      Histo->testinghitmap->Fill(Telescope[id]->Solution[9].Xpos, Telescope[id]->Solution[9].Ypos);
      if (id == 4) Histo->testWhit->Fill(Telescope[id]->Solution[9].ifront,Telescope[id]->Solution[9].iback);

    }

    //Front vs back and high vs low for gobbi
    if (Telescope[id]->Front.Nstore ==1 && Telescope[id]->Back.Nstore ==1)
    {
      float Elow_R = Telescope[id]->Front.Order[0].energylowR;
      int strip = Telescope[id]->Front.Order[0].strip;

      //Same thing here
      //Histo->GobbiFrontVsBack[id][strip]->SetPoint(Gfvb_cnt[id], Elow_R,Telescope[id]->Front.Order[0].energy);

      //Histo->GobbiBackHvsL[id][strip]->SetPoint(Gfvb_cnt[id],Telescope[id]->Back.Order[0].energylowR,Telescope[id]->Back.Order[0].energy); 

      Gfvb_cnt[id]++;
    }
    //continue;

    //Track num of coincidence events with just Si, CsI, both, and later matched solutions
    if (S800_results.trig_coin == true)
    {
      if (Telescope[id]->CsI.Nstore > 0 && Telescope[id]->Front.Nstore <= 0 && Telescope[id]->Back.Nstore <= 0) S800coinc_CsI = true;
      if (Telescope[id]->CsI.Nstore <= 0 && Telescope[id]->Front.Nstore > 0 && Telescope[id]->Back.Nstore > 0) S800coinc_Si = true;
      if (Telescope[id]->CsI.Nstore > 0 && Telescope[id]->Front.Nstore > 0 && Telescope[id]->Back.Nstore > 0) S800coinc_both = true;
    }

    //Place this after the alpha calibrations, f vs b, and high vs low
    //Since we have E-CsI only, no point in solutions without CsI. Continue if none
    //TODO Do we want to save Si solutions without CsI? No use for inv mass but could be useful for debugging
    if (Telescope[id]->CsI.Nstore < 1) continue;

		//Fill unmatched PSD spectra
		for (int i=0;i<Telescope[id]->CsI.Nstore;i++)
		{
      if (Telescope[id]->CsI.Order[i].qdcflag == 1)
      {
			  Histo->PSD_CsI_unmatched[id][i]->Fill(Telescope[id]->CsI.Order[i].energyR,Telescope[id]->CsI.Order[i].qdc);
      }		
    }

    //Make sure solutions are reset to 0 for each tele
    Telescope[id]->Nsolution = 0;
    Telescope[id]->CsINsolution = 0;

   //Want to fill f vs b unmatched strips here to look for things that will be matched out
   //No low gain strips until calibrations are done 
   if (Telescope[id]->Front.Nstore > 0 && Telescope[id]->Back.Nstore > 0)
   {
    float UMenergy;
    float UMbenergy;
    UMenergy = Telescope[id]->Front.Order[0].energy;
    UMbenergy = Telescope[id]->Back.Order[0].energy;
    //Histo->GFBE[id]->Fill(UMenergy,UMbenergy);
   }

    //E + CsI events
    //Work from least complex to most complex, least complex is all have 1 hit with CsI not matching
    //Most complex is same CsI quad with multiple silicons, needs zline matching
    //TODO make sure that all F, B, CsI cases are handled properly
    int CsIN;
    CsIN = Telescope[id]->CsI.Nstore;
    int FrontN = Telescope[id]->Front.Nstore;
    int BackN = Telescope[id]->Back.Nstore;


    //TODO what's going on here? Looks like a way of saving CsI before matching, similar code after the matching
    if (Telescope[id]->CsI.Nstore >= 1)
    {

      if (Telescope[id]->CsI.Nstore > 1)
      {
        for (int i=0;i<Telescope[id]->CsI.Nstore-1;i++)
        {
          if (Telescope[id]->CsI.Order[i].strip == Telescope[id]->CsI.Order[i+1].strip) NSameCsI++;
          if (Telescope[id]->CsI.Order[i].strip == Telescope[id]->CsI.Order[i+1].strip) break;
        }
      }


      //look for the simple case first
      //This is fine
      if (Telescope[id]->CsI.Nstore == 1 && Telescope[id]->Front.Nstore == 1 && Telescope[id]->Back.Nstore == 1)
      {
        Nmatch = Telescope[id]->simpleECsI();
        NsimpleECsI += Nmatch; //increase total count
        multiECsI += Nmatch;   //increase current count
      }

      //TODO This might be complex enough. Make sure that multiHitECsI have dEE matching. Won't work at first
      //FIXME DEE matching for multihits won't work until you have established zlines, turn off in telescope.cpp
      else if (Telescope[id]->Front.Nstore >=1 && Telescope[id]->Back.Nstore >=1)//then look at the multihit case
      {
        Nmatch = Telescope[id]->multiHitECsI();
        NmultiECsI += Nmatch;
        multiECsI += Nmatch;
      }

    }


    //Saves events with only CsI (no Si) into an seperate solutions array
    if ((Telescope[id]->Front.Nstore < 1 || Telescope[id]->Back.Nstore < 1) && CsIN >= 1)
    {
      for (int i=0;i<CsIN;i++)
      {
        Telescope[id]->CsISolution[Telescope[id]->CsINsolution].energy = Telescope[id]->CsI.Order[i].energy;
        Telescope[id]->CsISolution[Telescope[id]->CsINsolution].energyR = Telescope[id]->CsI.Order[i].energyR;
        Telescope[id]->CsINsolution++;
      }
    }
  }


  //plot E vs dE bananas and hitmap of paired E, CsI events, currently no hists for that
  //TODO all solutions will have CsI, don't need to check for that
  //Only high gain until you get good low cals
  //cout << "here1 " << endl;
  for (int id=0;id<4;id++) 
  {
    //Useful for tracking solutions across S800 coincidence events
    if (Telescope[id]->Nsolution > 0 && S800_results.trig_coin == true) S800_complete = true;

    for (int isol=0;isol<Telescope[id]->Nsolution; isol++)
    {
      //fill in hitmap of gobbi
      Telescope[id]->position(isol); //calculates x,y pos, and lab angle
      Histo->xyhitmap->Fill(Telescope[id]->Solution[isol].Xpos, Telescope[id]->Solution[isol].Ypos);
      //Theta phi hit map
      Histo->tphitmap->Fill(Telescope[id]->Solution[isol].theta*180./acos(-1.)*cos(Telescope[id]->Solution[isol].phi),Telescope[id]->Solution[isol].theta*180./acos(-1.)*sin(Telescope[id]->Solution[isol].phi));


      //Fill in E-CsI plots
      //TODO Might need to change from high to low gain, need good cals and know where to change
      float Ener = Telescope[id]->Solution[isol].energyR; // Want to use raw energy for CsI, draw zlines once
      float dEner;
      dEner = Telescope[id]->Solution[isol].denergy;
      int icsi = Telescope[id]->Solution[isol].iCsI;
      float theta = Telescope[id]->Solution[isol].theta;

      Histo->DEE_CsI[id][icsi]->Fill(Ener, dEner); //Not angle corrected

      int chan = id*4 + icsi;
      Histo->CsI_Energy_R[chan]->Fill(Ener);

      //Gates on S800 coincidences
      if (S800_results.trig_coin == true) Histo->DEE_CsI_S800Coinc[id][icsi]->Fill(Ener, dEner);
      if (S800_results.Zbeam == 14 && S800_results.Abeam == 23) Histo->DEE_CsI_S800Coinc_Si23[id][icsi]->Fill(Ener, dEner);

      Histo->DEE_CsI_0deg[id][icsi]->Fill(Ener, dEner*cos(theta)); //angle corrected
			Histo->PSD_CsI_matched[id][icsi]->Fill(Telescope[id]->Solution[isol].energyR,Telescope[id]->Solution[isol].qdc);

      //Fills raw CsI energy for center crystal hits
      // Chan     F     B     theta (deg)
      //  0      7,8   7,8        3.5
      //  1      7,8  23,24       7.0
      //  2     23,24 23,24       7.9
      //  3     23,24  7,8        5.2
      int fstrip = Telescope[id]->Solution[isol].ifront;
      int bstrip = Telescope[id]->Solution[isol].iback;
      if (fstrip == 7 || fstrip == 8 || fstrip == 23 || fstrip == 24)
      {
        if (bstrip == 7 || bstrip == 8 || bstrip == 23 || bstrip == 24)
        {
          Histo->CsI_Energy_R_center[chan]->Fill(Ener);
          //cout << " f b " << fstrip << " " << bstrip << " theta " << theta*180./acos(-1.) << endl;

        }
      }

      continue;
      //Fills hists for f vs b energy and matched strips
      Histo->GFBE_Matched[id]->Fill(Telescope[id]->Solution[isol].denergy,Telescope[id]->Solution[isol].benergy);

      //Looks at front and back strips vs CsI channels
      //Should be able to match telescope quadrants
      Histo->GFrontStrip_CsI[id]->Fill(icsi,Telescope[id]->Solution[isol].ifront);
      Histo->GBackStrip_CsI[id]->Fill(icsi,Telescope[id]->Solution[isol].iback);
     
      Histo->timediff_CsI[id][icsi]->Fill(Telescope[id]->Solution[isol].timediff);
    }
  }

  //calculate and determine particle identification PID in the Telescope
  int Pidmulti = 0;
  int nump_count = 0;
  flag1p = false;
  flag2p = false;
  flag3p = false;
  for (int id=0;id<4;id++) 
  {
    Pidmulti += Telescope[id]->getPID();

    //Tick up for protons to get num 2p
    for (int isol=0;isol<Telescope[id]->Nsolution; isol++)
    {

      if (Telescope[id]->Solution[isol].ipid == 1)
      {

        if (Telescope[id]->Solution[isol].iZ == 1 && Telescope[id]->Solution[isol].iA == 1)
        {
          num_protons++;          
          nump_count++;
        }
      }
    }
  }

  if (nump_count > 2) flag3p = true;

  if (nump_count == 1)
  {
    flag1p = true;
    num1p++;
  }
  if (nump_count == 2)
  {
    flag2p = true;
    num2p++;
  }

  //cout << num1p << endl;

  //calc sumEnergy,then account for Eloss in target, then set Ekin and momentum of solutions
  //Eloss files are loaded in Telescope
  //For hits in Gobbi, eloss through aluminum absorber and target
  //For S800 hit, eloss through BC400 and target

  for (int id=0;id<4;id++) 
  {
    Telescope[id]->calcEloss();

    for (int isol=0;isol<Telescope[id]->Nsolution;isol++)
    {
      //This is for histogramming the proton beam energy before target
      double pEkin = Telescope[id]->Solution[isol].Ekin;
      if (pEkin > 0) {
      Histo->CsI_Energy_protonBE[Telescope[id]->Solution[isol].iCsI + 4*id]->Fill(pEkin); }
    }

  }


  return multiECsI;
}


