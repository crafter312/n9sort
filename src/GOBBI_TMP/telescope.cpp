//Johnathan Phillips 2024_08_24
//class dealing with a single Si-CsI telescope
//       .-.      _______                             .  '  *   .  . '
//      {}``; |==|_______D                              .  * * -+-  
//      / ('        /|\                             .    * .    '  *
//  (  /  |        / | \                                * .  ' .  . 
//   \(_)_]]      /  |  \                            *   *  .   .
//                                                     '   *
#include "telescope.h"

//TODO
//    Erase all dE stuff, only have E-CsI
//    Implement simple and complex E-CsI matching
//    Erase all W silicon stuff
//    Add in loss through an aluminum absorber

using namespace std;

//constructor
telescope::telescope(bool S800)
{
 //Verify these numbers
  SiWidth = 6.42;
  SiFrame = 7.237;

  //switch which loss file is used depending on the target.

  if (S800 == false) Allosses = new CLosses(7,"_Al.loss"); //don't need Al loss for S800
  cout << "here2" << endl;

  Ran = CRandom::instance();
  Tthick = targthick::instance();

  cout << "here3" << endl;
  if (S800 == false) Targlosses = new CLosses(7,Tthick->TargType); //S800 targ loss in det.cpp
  cout << "Target Type: " << Tthick->TargType << endl;

  //CsI calibrations for clusters from proton calibrations
  string name = "cal/csi/deuteron.cal";
  calCsi_d = new calibrate(1,16,name,1,false);
  //cout << calCsi_d->Ntele << " " << calCsi_d->Nstrip << endl;
  name = "cal/csi/triton.cal";
  calCsi_t = new calibrate(1,16,name,1,false);
  //name = "Cal/csi/He3.cal";
  //calCsi_3He = new calibrate(1,16,name,1);
  name = "cal/csi/alpha.cal";
  calCsi_Alpha = new calibrate(1,16,name,1,false);

  
}

//destructor
telescope::~telescope()
{
  delete Targlosses;
  delete Allosses;
  delete Ran;
  
  for (int i=0;i<4;i++)  delete PidECsI[i];
  delete calCsi_d;
  delete calCsi_t;
  delete calCsi_Alpha;
}

//inialization
void telescope::init(int id0)
{
  id = id0;
  //-TODO check params for 35 mm Gobbi hole
  //     ____  Front view of Gobbi
  //    |    |____                  
  //    | 4  |    |              
  //   _|____| 1  |    __
  //--|--->|_|____|-->|5 | TODO could use telescope 5 for S800 solutions. If so, don't loop over it until loading sols
  //  | 3  |    |     |__|
  //  |____| 2  |
  //       |____|
  //TODO silicons may be slightly off centered, get actual positions. Also now at 35 mm hole, change centers
  //float const XcenterA[4] = {4.624,2.631,-4.624,-2.631};
  //float const YcenterA[4] = {2.631,-4.624,-2.631,4.624};
  float const XcenterA[4] = {1.5 + SiFrame/2.,SiFrame/2. - 1.5,-1.5 - SiFrame/2.,-SiFrame/2. + 1.5};
  float const YcenterA[4] = {SiFrame/2. - 1.5,-1.5 - SiFrame/2.,-SiFrame/2. + 1.5,1.5 + SiFrame/2.};
  Xcenter = XcenterA[id];
  Ycenter = YcenterA[id];
  

  ostringstream outstring;
  
  for (int i=0;i<4;i++)
  {
    outstring.str("");
    outstring << "pid_quad" << id+1 << "_CsI"<< i+1;
    PidECsI[i] = new pid(outstring.str(), false); //don't use S800 zline filepath
  }


}

void telescope::SetTarget(double dist, float thick)
{
  for (int i=0;i<10;i++){ Solution[i].SetTargetDistance(dist); }
  TargetThickness = thick;
}

void telescope::reset()
{
  multFront = 0;
  multBack = 0;
  multCsI = 0;

  Front.reset();
  tempFront.reset();

  Back.reset();
  tempBack.reset();

  CsI.reset();


  //possible to loop to Nsolution (but safer to just loop over all possible solution)
  for (int i=0; i<10; i++){ Solution[i].reset();}

  Nsolution = 0;

  //Reset CsI hits
  for (int i=0; i<10; i++){ CsISolution[i].reset();}

  CsINsolution = 0;
}

//????????
void telescope::Reduce()
{
  multFront = Front.Reduce("F");
  multBack = Back.Reduce("B");
}



//  .--.      .-'.      .--.      .--.      .--.      .--.      .`-.      .--.
//:::::.\::::::::.\::::::::.\::::::::.\::::::::.\::::::::.\::::::::.\::::::::.\
//'      `--'      `.-'      `--'      `--'      `--'      `-.'      `--'      `
//Here are 3 subroutines that are called to figure out how to go from detector 
//data to solutions that track a single particle.
//  -testingHitE()
//  -SimpleECsI()
//  -multiHitECsI() + loop()
//  Big one is multiHitdEECsI, acts as a hub to sort mixed multi events and send them to other subroutines
//  -multiHitdEECsI() -> sends events to all subroutines except testingHitE


// subroutine to identify a position map from alpha calibrations 
// it stores the answer in the last solution.
int telescope::testingHitE()
{
  Solution[9].energy = Front.Order[0].energy;
  Solution[9].energylow = Front.Order[0].energylow;
  Solution[9].energylowR = Front.Order[0].energylowR;
  Solution[9].energyR = Front.Order[0].energyR;
  Solution[9].benergy = Back.Order[0].energy;
  Solution[9].benergyR = Back.Order[0].energyR;
  Solution[9].denergy = -1;
  Solution[9].denergyR = -1;;
  Solution[9].ifront = Front.Order[0].strip;
  Solution[9].iback = Back.Order[0].strip;
  Solution[9].ide = -1;
  Solution[9].iCsI = -1;
  Solution[9].itele = id;
  Solution[9].timediff = -1000000.0;
  Solution[9].isSiCsI = false;

  Nsolution = 0; //don't determine there is any solutions.
  return 0;
}

// Simple E-CsI matching for 1 F, B, and CsI
int telescope::simpleECsI()
{
  //here we check the CsI is behind the x-y potions of the E
  //if it doesn't match, it is removed so there is a chance it is a dE-E event

  if (Front.Order[0].strip <= 16 && Back.Order[0].strip <= 16 && 
      CsI.Order[0].strip !=0)
  {
    Nsolution = 0;

    //Store CsI (no Si) into an array for later use. Then remove
    CsISolution[CsINsolution].energy = CsI.Order[0].energy;
    CsISolution[CsINsolution].energyR = CsI.Order[0].energyR;
    CsINsolution++;
    CsI.Remove(0);
    return 0;
  }
  else if (Front.Order[0].strip <= 16 && Back.Order[0].strip >= 15 && 
      CsI.Order[0].strip !=1)
  {
    Nsolution = 0;

    CsISolution[CsINsolution].energy = CsI.Order[0].energy;
    CsISolution[CsINsolution].energyR = CsI.Order[0].energyR;
    CsINsolution++;
    CsI.Remove(0);
    return 0;
  }
  else if (Front.Order[0].strip >= 15 && Back.Order[0].strip >= 15 && 
      CsI.Order[0].strip !=2)
  {
    Nsolution = 0;

    CsISolution[CsINsolution].energy = CsI.Order[0].energy;
    CsISolution[CsINsolution].energyR = CsI.Order[0].energyR;
    CsINsolution++;
    CsI.Remove(0);
    return 0;
  }
  else if (Front.Order[0].strip >= 15 && Back.Order[0].strip <= 16 && 
      CsI.Order[0].strip !=3)
  {
    Nsolution = 0;

    CsISolution[CsINsolution].energy = CsI.Order[0].energy;
    CsISolution[CsINsolution].energyR = CsI.Order[0].energyR;
    CsINsolution++;
    CsI.Remove(0);
    return 0;
  }

  //this won't work until your calibrations are good. but it should be turned on.
  if (fabs(Front.Order[0].energy - Back.Order[0].energy) > 10.) 
  {
    Nsolution = 0;

    CsISolution[CsINsolution].energy = CsI.Order[0].energy;
    CsISolution[CsINsolution].energyR = CsI.Order[0].energyR;
    CsINsolution++;
    CsI.Remove(0);
    return 0;
  }

  float timediff = CsI.Order[0].time - Front.Order[0].time;
  Solution[Nsolution].energy = CsI.Order[0].energy;
  Solution[Nsolution].energyR = CsI.Order[0].energyR;
  Solution[Nsolution].energylow = CsI.Order[0].energy; //The later code uses high gain and low gain for gobbi
  Solution[Nsolution].energylowR = CsI.Order[0].energyR; //To make the code simpler, just clone CsI high into low gain
  Solution[Nsolution].benergy = Back.Order[0].energy;
  Solution[Nsolution].benergyR = Back.Order[0].energyR;
  Solution[Nsolution].denergy = Front.Order[0].energy;
  Solution[Nsolution].denergylow = Front.Order[0].energylow;
  Solution[Nsolution].denergyR = Front.Order[0].energyR;
  Solution[Nsolution].qdc = CsI.Order[0].qdc;

  Solution[Nsolution].ifront = Front.Order[0].strip;
  Solution[Nsolution].iback = Back.Order[0].strip;
  Solution[Nsolution].ide = -1;
  Solution[Nsolution].iCsI= CsI.Order[0].strip;
  Solution[Nsolution].itele = id; 
  Solution[Nsolution].CsITime = CsI.Order[0].time;
  Solution[Nsolution].timediff = timediff;
  Solution[Nsolution].fTime = Front.Order[0].time;
  Solution[Nsolution].isSiCsI = true;
  Solution[Nsolution].itele = id;
  //cout << "telescope::simple iCsI=" <<  Solution[0].iCsI << " itele = " << id << "solution = 0" << endl;
  Nsolution += 1;
  
  return 1;
}

//****************************************************
//recursive subroutine  used for multihit subroutine
//this is a complicated subroutine so i have a mega comment labled multihit_comment.h 
//stored in this directory. It will step you through a multihit Si-Si event
//TODO Ne16 version also matched dE events, only care about F and B now. Erase dE and make sure it works
void telescope::loop(int depth)
{
  //only work in this section if we are at the max level of recursion
  if (depth == NestDim) 
  {
    //as an example for NestDim = 2
    //for first time through on loop(2) nestarray = {0,1} so we check
    //if highest energy delta matches with highest energy Front. Then 
    //check if second highest delta matches highest energy Front.

    float de = 0.;
    for (int i=0;i<NestDim;i++)
    {

      //difference in Front and back energy is how to match FrontE and BackE - need two types for Gobbi
      //TODO this needs high and low calibrations to work
      //For E < 30
      if (Front.Order[i].energy < 30.) de += abs(Back.Order[NestArray[i]].energy-Front.Order[i].energy);
      //For E >= 30
      if (Front.Order[i].energy >= 30.) de += abs(Back.Order[NestArray[i]].energylow-Front.Order[i].energylow);
    }

    //here if it is the lowest total difference in strip# or energy, it is saved in
    //arrayB (for matching front-back)
    if (de < deMin)
    {
      deMin = de;
      for (int i=0;i<NestDim;i++) {arrayB[i] = NestArray[i];}
    }
    return;
  }

  //this section handles how deep we go into the recursion loop.
  //the key to this section if figuring out what NestArrays to check.
  //for NestDim=2, we want to check {0,1} and {1,0}
  for (int i=0;i<NestDim;i++)
  {
    NestArray[depth] = i;
    int leave = 0;
    //when matching we want to skip items already matched 
    for (int j=0;j<depth;j++)
    {
      if (NestArray[j] == i)
      {
        leave =1;
        break; 
      }
    }
    if (leave) continue;
    loop(depth+1);
  }
}


//TODO this was all new code and needs to be checked over thoroughly
//***************************************************
//modification of multiHitdEE() to match E-CsI events
int telescope::multiHitECsI()
{
  int Ntries = min(Front.Nstore,Back.Nstore);
  //Ntries = min(Ntries,CsI.Nstore);

  if (Ntries > 4)
    Ntries = 4;
  if (Ntries <= 0)
    return 0;

  NSisolution = 0;
  for (NestDim = Ntries;NestDim>0;NestDim--)
  {
    dstripMin = 1000;
    deMin = 10000.;
    
    //look for best Front/Back matching
    //solutions are stored in arrayB
    loop(0);
    
    //check to see if best possible solution is reasonable
    int leave = 0;
    for (int i=0;i<NestDim;i++)
    {
      if (fabs(Back.Order[arrayB[i]].energy - Front.Order[i].energy) > 10.) 
      {
        leave = 1;
        break;
      }
    }
    if (leave) continue;
    NSisolution = NestDim;
  }
  //save some time and just exit if there are no Si Front/Back solutions
  if (NSisolution == 0) return 0;

  //now assign each of these solutions a Csi detector location 
  int mult[4]={0};  //array for multipility of Si solution for each Csi
  int sil[4][10];   //contains a lits of solutions for each Csi

  //look at all the Front/Back solutions and see how many are on each CsI
  for (int i=0;i<NSisolution;i++)
  {
    int ifront = Front.Order[i].strip;
    int iback = Back.Order[arrayB[i]].strip;
    for (int icsi=0;icsi<4;icsi++)
    {
      if (ifront >= FrontLow[icsi] &&
          ifront <= FrontHigh[icsi] &&
          iback  >= BackLow[icsi]  &&
          iback  <= BackHigh[icsi])
      {
        sil[icsi][mult[icsi]] = i;
        mult[icsi]++;
        break;
      }
    }
  }

  //make array of detect csi energies
  float energy[4]={-1.};
  float energyR[4]={-1.};
  short order[4]={-1};

  //store the CsI raw energy info in an array that corresponds to the position it is in
  for (int i=0;i<CsI.Nstore;i++)
  {
    energy[CsI.Order[i].strip] = CsI.Order[i].energy;
    order[CsI.Order[i].strip] = i;
  }

  //Define the remove constants as some index that won't happen
  int removeFirst = -1;
  int removeSecond = -1;

  //loop over csi location
  for (int icsi = 0;icsi<4;icsi++)
  {
    //no solution for this location, ignore
    if (mult[icsi] == 0) continue;

    //FIXME DEE matching won't work until you have good zlines. Turn off at start
    //If more than 1 si solution for a single CsI, check if it falls in a zline
    //CsI can only fire once within readout
    //Needed for events with mixed dE and CsI in the same quad
    //Can only accept one solution, don't allow it to accept both
    else if (mult[icsi] > 1)
    {
      continue;
      /*for (int i=0;i<mult[icsi];i++)
      {
        int ii = sil[icsi][i];
        //Do zline check
        int zCheck = 0;

        zCheck = PidECsI[id]->getEGate(CsI.Order[0].energyR,Front.Order[ii].energy);


        if (zCheck == 0) continue;
        else 
        { //need to fill stuff here using the correct "ii" index
          Solution[Nsolution].energy = energy[icsi];
          Solution[Nsolution].energyR = CsI.Order[order[icsi]].energyR;
          Solution[Nsolution].energylow = energy[icsi];
          Solution[Nsolution].energylowR = CsI.Order[order[icsi]].energyR;
          Solution[Nsolution].denergy = Front.Order[ii].energy;
          Solution[Nsolution].denergylow = Front.Order[ii].energylow;
          Solution[Nsolution].denergyR = Front.Order[ii].energyR;
          Solution[Nsolution].benergy = Back.Order[arrayB[ii]].energy;
          Solution[Nsolution].benergylow = Back.Order[arrayB[ii]].energylow;
          Solution[Nsolution].benergyR = Back.Order[arrayB[ii]].energyR;
      		Solution[Nsolution].qdc = CsI.Order[order[icsi]].qdc;

          Solution[Nsolution].ifront = Front.Order[ii].strip;
          Solution[Nsolution].iback = Back.Order[arrayB[ii]].strip;
          Solution[Nsolution].ide = -1;
          Solution[Nsolution].iCsI = icsi;
          Solution[Nsolution].itele = id;
          Solution[Nsolution].isSiCsI = true;
          Solution[Nsolution].CsITime = CsI.Order[order[icsi]].time;
          float timediff = CsI.Order[order[icsi]].time - Front.Order[ii].time;
          Solution[Nsolution].timediff = timediff;
          Solution[Nsolution].fTime = Front.Order[ii].time;
          Nsolution++;

          Front.Order[ii].CsIFlag = 1;
          Back.Order[arrayB[ii]].CsIFlag = 1;

          break; //break out of loop, accept only one solution per quad. Allows for small chance of losing dE solution
          //I could find a more elegant solution using strip matching for dE and base it on a best score
        }
      }*/
    }

    //more than one si solution for a single Csi location ignore
    //CsI can only fire once within readout
    //else if (mult[icsi] > 1) continue;
    // CsI energy < 0 should not happen
    //but just in case ignore
    else if(energy[icsi] <= 0.) continue;
    
    else
    {

      int ii = sil[icsi][0];
      Solution[Nsolution].energy = energy[icsi];
      Solution[Nsolution].energyR = CsI.Order[order[icsi]].energyR;
      Solution[Nsolution].energylow = energy[icsi];
      Solution[Nsolution].energylowR = CsI.Order[order[icsi]].energyR;
      Solution[Nsolution].denergy = Front.Order[ii].energy;
      Solution[Nsolution].denergylow = Front.Order[ii].energylow;
      Solution[Nsolution].denergyR = Front.Order[ii].energyR;
      Solution[Nsolution].benergy = Back.Order[arrayB[ii]].energy;
      Solution[Nsolution].benergylow = Back.Order[arrayB[ii]].energylow;
      Solution[Nsolution].benergyR = Back.Order[arrayB[ii]].energyR;
      Solution[Nsolution].qdc = CsI.Order[order[icsi]].qdc;

      Solution[Nsolution].ifront = Front.Order[ii].strip;
      Solution[Nsolution].iback = Back.Order[arrayB[ii]].strip;
      Solution[Nsolution].ide = -1;
      Solution[Nsolution].iCsI = icsi;
      Solution[Nsolution].itele = id;
      Solution[Nsolution].isSiCsI = true;
      Solution[Nsolution].CsITime = CsI.Order[order[icsi]].time;
      float timediff = CsI.Order[order[icsi]].time - Front.Order[ii].time;
      Solution[Nsolution].timediff = timediff;
      Solution[Nsolution].fTime = Front.Order[ii].time;
      Nsolution++;

      Front.Order[ii].CsIFlag = 1;
      Back.Order[arrayB[ii]].CsIFlag = 1;
      
    }
  }

  return Nsolution;
}

//*************************************
//finds particle identification - checks to see if particle is inside of z - bananas  
int telescope::getPID()
{
  int pidmulti = 0;
  for (int isol=0; isol<Nsolution; isol++)
  {

    Solution[isol].ipid = 0;
    bool isSiCsI = Solution[isol].isSiCsI;

    //Finds angle corrected angle to get PID from DEE zlines
    float energy;
    float energyR;
    float denergy;
    float angle;
    angle = Solution[isol].theta;
    //Gobbi has a mix of high and low gain for E and dE, telescope 1 has E back energy.
    //Angle correct dE at the end to preserve high and low gain dE without too much code
    //Use energyR for CsI silicon dE, high gain, no angle correction
    //TODO need good high and low cals
    //if (denergy < 25.) denergy = Solution[isol].denergy;
    //if (denergy >= 25.) denergy = Solution[isol].denergylow;
    denergy = Solution[isol].denergy;

    denergy = denergy*cos(angle); //angle correct dE now
    energyR = Solution[isol].energyR;
    bool FoundPid = false;
    if (isSiCsI)
    {
      //use raw energy for CsI PID
      //cout << Solution[isol].iCsI << endl;
      FoundPid = PidECsI[Solution[isol].iCsI]->getPID(energyR, denergy);
      //Pid->Z = PidECsI[Solution[isol].iCsI]->Z;
      //Pid->A = PidECsI[Solution[isol].iCsI]->A;
      //Pid->mass = PidECsI[Solution[isol].iCsI]->mass;
      //if(Pid->Z == 1 && Pid->A == 1) cout << "proton Si-CsI" << endl;

    }
    else
    {
      cout << "No CsI solution, shouldn't be in this subroutine " << endl; 
    }

    //no particle id is found
    if (!FoundPid) continue;
    else pidmulti++;
    if (pidmulti == 2) cout << pidmulti << endl;
    Solution[isol].ipid = 1; //this can be adapted to be different values later
    Solution[isol].iZ = PidECsI[Solution[isol].iCsI]->Z;
    Solution[isol].iA = PidECsI[Solution[isol].iCsI]->A;
    

    Solution[isol].mass = PidECsI[Solution[isol].iCsI]->getMass(PidECsI[Solution[isol].iCsI]->Z,PidECsI[Solution[isol].iCsI]->A); //we want mass in energy units not AMU
    //cout << Solution[isol].mass << endl;

    //take proton equivalent energies to light equivalent
    if (isSiCsI)
    {
      int csid = Solution[isol].iCsI + 4*id;
      //Solution[isol].energy = light2energy(Pid->Z, Pid->A, csid, Solution[isol].energy); turned off light to energy
    }

  }
  return pidmulti;
}
//********************************************************************

int telescope::calcEloss()
{
  for (int isol=0; isol<Nsolution; isol++)
  {
    //need PID to calculate energy loss
    if (!Solution[isol].ipid)
    {
      Solution[isol].Ekin = 0;
      return 0;
    }

    //kinetics calc, add Delta and energy for total energy
    float sumEnergy = Solution[isol].denergy + Solution[isol].energy;

//*****************************************
    float pc_before = sqrt(pow(sumEnergy+Solution[isol].mass,2) - pow(Solution[isol].mass,2));
    float velocity_before = pc_before/(sumEnergy+Solution[isol].mass);

    //  HEY YOU!    Make sure that only Gobbi events go here. Gobbi events have loss through the target and an Al absorber
    
    //Three target thickness for 9Be

    //Losses through 5 mm Al absorber. Only for Gobbi hits
    float Althick = 1351.; //mg/cm^2
    Althick = Althick/cos(Solution[isol].theta);
    float Alein = 0;
    Alein = Allosses->getEin(sumEnergy,Althick,Solution[isol].iZ,Solution[isol].mass/m0);
    //cout << "AL " << sumEnergy << " " << Alein << endl;
    //Losses through 9Be target
    float thick = TargetThickness/2./cos(Solution[isol].theta);
    //float thick = 369.6/cos(Solution[isol].theta); //For the Be proton cal run (2mm targ)
    //cout << "thick 0 deg full " << TargetThickness << " half " << thick << endl;
    //Targlosses = new CLosses(8,TargType); //have to initialize here to include different targets
    float ein = Targlosses->getEin(Alein,thick,Solution[isol].iZ,Solution[isol].mass/m0);
    //cout << "Be " << Alein << " " << ein << endl;
    //ein = ein*0.99; //REDUCE THE CALIBRATION BY SOME FACTOR
    Solution[isol].Ekin = ein;

    //calc momentum vector, energyTot, and velocity
    Solution[isol].getMomentum();
  }

  return 1;
}


//*******************************************************************************
  //calculates the x-y position and angles in the array in cm
void telescope::position(int isol)
{
  float Xpos,Ypos;

  if (id == 0) 
  {
    Xpos = Xcenter + 
          (((double)Solution[isol].iback+Ran->Rndm())/32.-0.5)*SiWidth;
    Ypos = Ycenter +
          (((double)Solution[isol].ifront+Ran->Rndm())/32.-0.5)*SiWidth;
  }
  else if (id == 1)
  {
    Xpos = Xcenter + 
          (((double)Solution[isol].ifront+Ran->Rndm())/32.-0.5)*SiWidth;
    Ypos = Ycenter +
          (0.5-((double)Solution[isol].iback+Ran->Rndm())/32.)*SiWidth;
  }
  else if (id == 2)
  {
    Xpos = Xcenter + 
          (0.5-((double)Solution[isol].iback+Ran->Rndm())/32.)*SiWidth;
    Ypos = Ycenter +
          (0.5-((double)Solution[isol].ifront+Ran->Rndm())/32.)*SiWidth;
  }
  else if (id == 3)
  {
    Xpos = Xcenter + 
          (0.5-((double)Solution[isol].ifront+Ran->Rndm())/32.)*SiWidth;
    Ypos = Ycenter +
          (((double)Solution[isol].iback+Ran->Rndm())/32.-0.5)*SiWidth;
  }

  Solution[isol].Xpos = Xpos;
  Solution[isol].Ypos = Ypos;
  float theta = Solution[isol].angle();
}
//***********************************************************************

/* only need positionC if you don't want to randomly distribute event in a strip
//*******************************************************************************
  //calculates the x-y position in the array in cm
void telescope::positionC(int isol)
{
  float Xpos,Ypos;

  if (id == 0) 
  {
    Xpos = Xcenter +
         (((double)Solution[isol].iback+.5)/32.-0.5)*SiWidth;
    Ypos = Ycenter +
         (((double)Solution[isol].ifront+.5)/32.-0.5)*SiWidth;
  }
  else if (id == 1)
  {
    Xpos = Xcenter + 
          (((double)Solution[isol].ifront+.5)/32.-0.5)*SiWidth;
    Ypos = Ycenter +
          (0.5-((double)Solution[isol].iback+.5)/32.)*SiWidth;
  }
  else if (id == 2)
  {
    Xpos = Xcenter +
          (0.5-((double)Solution[isol].iback+.5)/32.)*SiWidth;
    Ypos = Ycenter +
          (0.5-((double)Solution[isol].ifront+.5)/32.)*SiWidth;
  }
  else if (id == 3)
  {
    Xpos = Xcenter +
          (0.5-((double)Solution[isol].ifront+.5)/32.)*SiWidth;
    Ypos = Ycenter +
          (((double)Solution[isol].iback+.5)/32.-0.5)*SiWidth;
  }
  Solution[isol].Xpos = Xpos;light2
  Solution[isol].Ypos = Ypos;
  float theta = Solution[isol].angle();
}
*/

//***********************************************************
// load front and back strip delimeters for each CsI crystal.
void telescope::load(int F0low, int F1low,int F2low, int F3low,
                 int F0hi,  int F1hi, int F2hi,  int F3hi,
                 int B0low, int B1low,int B2low, int B3low,
                 int B0hi,  int B1hi, int B2hi,  int B3hi)
{
  FrontLow[0] = F0low;
  FrontLow[1] = F1low;
  FrontLow[2] = F2low;
  FrontLow[3] = F3low;

  FrontHigh[0] = F0hi;
  FrontHigh[1] = F1hi;
  FrontHigh[2] = F2hi;
  FrontHigh[3] = F3hi;

  BackLow[0] = B0low;
  BackLow[1] = B1low;
  BackLow[2] = B2low;
  BackLow[3] = B3low;

  BackHigh[0] = B0hi;
  BackHigh[1] = B1hi;
  BackHigh[2] = B2hi;
  BackHigh[3] = B3hi;

}

//***************************************************
// converstion of equilivant proton energy to energy for a given isotope
//i.e. Z and A dependence of CsI light output
float telescope::light2energy(int Z, int A, int CsIhit, float energy)
{
  //double iniE = energy;
  //cout << endl << energy << endl;
  if(Z ==1)
  { 
    if(A ==2)
      {

      energy = calCsi_d->getEnergy(0,CsIhit,energy);
       }
    if(A ==3)
      energy = calCsi_t->getEnergy(0,CsIhit,energy);
  }
  else if(Z == 2)
  {
    if(A ==3)
      energy = calCsi_Alpha->getEnergy(0,CsIhit,energy);
    else if(A ==4)
      energy = calCsi_Alpha->getEnergy(0,CsIhit,energy);
    else 
    {
      //cout << "found no calib for " << Z << " " << A << endl; 
      return -1.;
    }
  }
  else 
  {
    cout << "found no calib for Z= " << Z << " A= " << A << " id =  " << id << " Csi = " << CsIhit << endl; 
    abort();
  }
  //if (iniE != energy) cout << Z << " " << A << " " << energy << endl;
  return energy;
}
//***********************************************************************


