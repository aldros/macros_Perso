#include <iomanip>
#include <iostream>
#include "TGeoManager.h"
#include "TMath.h"

// Create Matrix Unity
TGeoRotation *fGlobalRot = new TGeoRotation();

// Create a null translation
TGeoTranslation *fGlobalTrans = new TGeoTranslation();
TGeoRotation *fRefRot = NULL;

TGeoManager*  gGeoMan = NULL;

Double_t fThetaX = 0.;
Double_t fThetaY = 0.;
Double_t fThetaZ = 0.;
Double_t fPhi   = 0.;
Double_t fTheta = 0.;
Double_t fPsi   = 0.;
Double_t fX = 0.;
Double_t fY = 0.;
Double_t fZ = 0.;
Bool_t fLocalTrans = kFALSE;

TGeoCombiTrans* GetGlobalPosition(TGeoCombiTrans *fRef);

void create_foot_geo(const char* geoTag="2022")
{
  // --------------------------------------------------------------------------
  // Configurable geometry for the FOOT Detectors.
  // Use this macro to create root files with the different configurations 
  // and positions/angles of the silicon detectors.
  //
  // Execute macro:  root -l
  //                 .L create_foot_2022_geo.C
  //                 create_foot_geo()
  // --------------------------------------------------------------------------

  fGlobalTrans->SetTranslation(0.0,0.0,0.0);

  // -------   Load media from media file   -----------------------------------
  FairGeoLoader*    geoLoad = new FairGeoLoader("TGeo","FairGeoLoader");
  FairGeoInterface* geoFace = geoLoad->getGeoInterface();
  TString geoPath = gSystem->Getenv("VMCWORKDIR");
  TString medFile = geoPath + "/geometry/media_r3b.geo";
  geoFace->setMediaFile(medFile);
  geoFace->readMedia();
  gGeoMan = gGeoManager;
  // --------------------------------------------------------------------------



  // -------   Geometry file name (output)   ----------------------------------
  TString geoFileName = geoPath + "/geometry/foot_";
  geoFileName = geoFileName + geoTag + ".geo.root";
  // --------------------------------------------------------------------------

  // -----------------   Get and create the required media    -----------------
  FairGeoMedia*   geoMedia = geoFace->getMedia();
  FairGeoBuilder* geoBuild = geoLoad->getGeoBuilder();

  FairGeoMedium* mAir      = geoMedia->getMedium("Air");
  if ( ! mAir ) Fatal("Main", "FairMedium Air not found");
  geoBuild->createMedium(mAir);
  TGeoMedium* pMedAir = gGeoMan->GetMedium("Air");
  if ( ! pMedAir ) Fatal("Main", "Medium Air not found");

  FairGeoMedium* mVac      = geoMedia->getMedium("vacuum");
  if ( ! mVac ) Fatal("Main", "FairMedium vacuum not found");
  geoBuild->createMedium(mVac);
  TGeoMedium* pMedVac = gGeoMan->GetMedium("vacuum");
  if ( ! pMedVac ) Fatal("Main", "Medium vacuum not found");
  
  FairGeoMedium* mSi      = geoMedia->getMedium("silicon");
  if ( ! mSi ) Fatal("Main", "FairMedium silicon not found");
  geoBuild->createMedium(mSi);
  TGeoMedium* pMedSi = gGeoMan->GetMedium("silicon");
  if ( ! pMedSi ) Fatal("Main", "Medium silicon not found");

  FairGeoMedium* mCopper      = geoMedia->getMedium("copper");
  if ( ! mCopper ) Fatal("Main", "FairMedium copper not found");
  geoBuild->createMedium(mCopper);
  TGeoMedium* pMedCu = gGeoMan->GetMedium("copper");
  if ( ! pMedCu ) Fatal("Main", "Medium copper not found");
  
  FairGeoMedium* mAl      = geoMedia->getMedium("aluminium");
  if ( ! mAl ) Fatal("Main", "FairMedium aluminium not found");
  geoBuild->createMedium(mAl);
  TGeoMedium* pMedAl = gGeoMan->GetMedium("aluminium");
  if ( ! pMedAl ) Fatal("Main", "Medium aluminium not found");

  FairGeoMedium* mFe      = geoMedia->getMedium("iron");
  if ( ! mFe ) Fatal("Main", "FairMedium iron not found");
  geoBuild->createMedium(mFe);
  TGeoMedium* pMedFe = gGeoMan->GetMedium("iron");
  if ( ! pMedFe ) Fatal("Main", "Medium iron not found");
  
  FairGeoMedium* mPolyethylene      = geoMedia->getMedium("polyethylene");
  if ( ! mPolyethylene ) Fatal("Main", "FairMedium polyethylene not found");
  geoBuild->createMedium(mPolyethylene);
  TGeoMedium* pMedPolyethylene = gGeoMan->GetMedium("polyethylene");
  if ( ! pMedPolyethylene ) Fatal("Main", "Medium polyethylene not found");
  // --------------------------------------------------------------------------



  // --------------   Create geometry and top volume  -------------------------
  gGeoMan = (TGeoManager*)gROOT->FindObject("FAIRGeom");
  gGeoMan->SetName("TRAgeom");
  TGeoVolume* top = new TGeoVolumeAssembly("TOP");
  gGeoMan->SetTopVolume(top);
  // --------------------------------------------------------------------------



  // out-of-file geometry definition
  Double_t dx,dy,dz;
  Double_t rmin, rmax, rmin1, rmax1, rmin2, rmax2;
  Double_t thx, phx, thy, phy, thz, phz;
  Double_t  phi1, phi2;
  Double_t tx,ty,tz;
  
  TGeoRotation * zeroRot = new TGeoRotation; //zero rotation
  TGeoCombiTrans * tZero = new TGeoCombiTrans("tZero", 0., 0., 0., zeroRot);
  tZero->RegisterYourself();
  
  
  //-------------------------------------------------------------------
  
  // Fill Chamber: Vacuum or Air. Needed still: an external call interface for choosing which.
  TGeoMedium * pMedFill=pMedVac;
  //pMedFill = new TGeoMedium("Fill_Air", numed,pMat2, par);
  //pMedFill = (TGeoMedium*) pMedAir->Clone();
  //pMedFill->SetName("Fill_Air");
  //  pMedFill = (TGeoMedium*) pMedVac->Clone();
  //  pMedFill->SetName("Fill_Vacuum");
  
  //-------------------------------------------------------------------
  
  
  // Shape: World type: TGeoBBox
  TGeoVolume* pWorld = gGeoManager->GetTopVolume();
  pWorld->SetVisLeaves(kTRUE);

  //------------------ PCBs------------------------------------------
  TGeoShape* padle_h_box1 = new TGeoBBox("padle_h_box1",
					 7.25, 
					 9.1, 
					 0.0150);
  TGeoShape* padle_h_box2 = new TGeoBBox("padle_h_box2",
					 5., 
					 5., 
					 0.0150);
  TGeoShape* padle_h_box3 = new TGeoBBox("padle_h_box3",
					 9.1, 
					 7.25, 
					 0.0150);
  
  //Shift the Si position
  TGeoCombiTrans * tShift = new TGeoCombiTrans("tShift", 0.0, -1.48, 0.015, zeroRot);
  tShift->RegisterYourself();

  TGeoCombiTrans * tShift2 = new TGeoCombiTrans("tShift2", 0.0, 1.48, 0.015, zeroRot);
  tShift2->RegisterYourself();

  TGeoCombiTrans * tShift3 = new TGeoCombiTrans("tShift3", -1.48, 0.0, 0.015, zeroRot);
  tShift3->RegisterYourself();

  TGeoCombiTrans * tShift4 = new TGeoCombiTrans("tShift4", 1.48, 0.0, 0.015, zeroRot);
  tShift4->RegisterYourself();
  
  // Create a composite shape
  TGeoCompositeShape *PCB = new TGeoCompositeShape("diffbox", "padle_h_box1:tShift - padle_h_box2");
  TGeoCompositeShape *PCB2 = new TGeoCompositeShape("diffbox2", "padle_h_box1:tShift2 - padle_h_box2");
  TGeoCompositeShape *PCB3 = new TGeoCompositeShape("diffbox3", "padle_h_box3:tShift3 - padle_h_box2");
  TGeoCompositeShape *PCB4 = new TGeoCompositeShape("diffbox4", "padle_h_box3:tShift4 - padle_h_box2");
  TGeoVolume *PCBvol1 = new TGeoVolume("PCB", PCB, pMedPolyethylene);
  PCBvol1->SetVisLeaves(kTRUE);
  PCBvol1->SetLineColor(kGreen+2);
  TGeoVolume *PCBvol2 = new TGeoVolume("PCB2", PCB2, pMedPolyethylene);
  PCBvol2->SetVisLeaves(kTRUE);
  PCBvol2->SetLineColor(kGreen+2);
  TGeoVolume *PCBvol3 = new TGeoVolume("PCB3", PCB3, pMedPolyethylene);
  PCBvol3->SetVisLeaves(kTRUE);
  PCBvol3->SetLineColor(kGreen+2);
  TGeoVolume *PCBvol4 = new TGeoVolume("PCB4", PCB4, pMedPolyethylene);
  PCBvol4->SetVisLeaves(kTRUE);
  PCBvol4->SetLineColor(kGreen+2);
  
  //**************************************************************//
  //*********************   Si Sensors ******************************//
  //***************************************************************/
  
  // Si Shape & volume: TraBox type: TGeoBBox
  dx = 5.; //0.5*14.50000; //5.
  dy = 5.; //0.5*(10.00000+2.*2.62); //5.
  dz = 0.0075000;
  // Volume: TraLog
  TGeoVolume *TraLog = gGeoManager->MakeBox("TraLog",pMedSi,dx,dy,dz);
  TraLog->SetVisLeaves(kTRUE);
  TraLog->SetLineColor(33);
  
  //TRANSFORMATION MATRICES
  dx = -6.488000;
  dy = 0.000000;
  dz = sqrt(8.47*8.47-dx*dx);
  // Rotation:
  thx = 33.400000;    phx = 0.000000;
  thy = 90.000000;    phy = 90.000000;
  thz = 0.000000;    phz = 0.000000;
  //TGeoRotation *      pMatrix3 = new TGeoRotation("",thx,phx,thy,phy,thz,phz);
  TGeoRotation *      pMatrix3 = new TGeoRotation();
  pMatrix3->SetAngles(90,-56.6,0);
  TGeoCombiTrans*  pMatrix2 = new TGeoCombiTrans("", dx,dy,dz,pMatrix3);

  //TRANSFORMATION MATRICES
  dx = -6.488000;
  dy = 0.000000;
  dz = sqrt(8.47*8.47-dx*dx)+0.1;
  // Rotation:
  thx = 33.400000;    phx = 0.000000;
  thy = 90.000000;    phy = 90.000000;
  thz = 0.000000;    phz = 0.000000;
  //TGeoRotation *      pMatrix3 = new TGeoRotation("",thx,phx,thy,phy,thz,phz);
  TGeoRotation *      pMatrix29 = new TGeoRotation();
  pMatrix29->SetAngles(90,-56.6,0);
  TGeoCombiTrans*  pMatrix28 = new TGeoCombiTrans("", dx,dy,dz,pMatrix29); 
  
  //Combi transformation:
  dx = -7.566000;
  dy = 0.000000;
  dz = sqrt(10.7*10.7-dx*dx);
  // Rotation:
  thx = 41.000000;    phx = 0.000000;
  thy = 90.000000;    phy = 90.000000;
  thz = 0.000000;    phz = 0.000000;
  //TGeoRotation *       pMatrix17 = new TGeoRotation("",thx,phx,thy,phy,thz,phz);
  TGeoRotation *       pMatrix17 = new TGeoRotation();
  pMatrix17->SetAngles(90,-49,0);
  TGeoCombiTrans*   pMatrix16 = new TGeoCombiTrans("", dx,dy,dz,pMatrix17);

  //Combi transformation:
  dx = -7.566000;
  dy = 0.000000;
  dz = sqrt(10.7*10.7-dx*dx)+0.1;
  // Rotation:
  thx = 41.000000;    phx = 0.000000;
  thy = 90.000000;    phy = 90.000000;
  thz = 0.000000;    phz = 0.000000;
  //TGeoRotation *       pMatrix17 = new TGeoRotation("",thx,phx,thy,phy,thz,phz);
  TGeoRotation *       pMatrix25 = new TGeoRotation();
  pMatrix25->SetAngles(90,-49,0);
  TGeoCombiTrans*   pMatrix24 = new TGeoCombiTrans("", dx,dy,dz,pMatrix25); 
	
  //Combi transformation:
  dx = 6.488000;
  dy = 0.000000;
  dz = sqrt(8.47*8.47-dx*dx);
  // Rotation:
  thx = -33.400000;   phx = 0.000000;
  thy = 90.000000;    phy = 90.000000;
  thz = 0.000000;    phz = 0.000000;
  //TGeoRotation *       pMatrix5 = new TGeoRotation("",thx,phx,thy,phy,thz,phz);
  TGeoRotation *       pMatrix5 = new TGeoRotation();
  pMatrix5->SetAngles(90,56.6,0); //56.6
  TGeoCombiTrans*   pMatrix4 = new TGeoCombiTrans("", dx,dy,dz,pMatrix5);

  //Combi transformation:
  dx = 6.488000;
  dy = 0.000000;
  dz = sqrt(8.47*8.47-dx*dx)+0.1;
  // Rotation:
  thx = -33.400000;    phx = 0.000000;
  thy = 90.000000;    phy = 90.000000;
  thz = 0.000000;    phz = 0.000000;
  //TGeoRotation *       pMatrix5 = new TGeoRotation("",thx,phx,thy,phy,thz,phz);
  TGeoRotation *       pMatrix27 = new TGeoRotation();
  pMatrix27->SetAngles(90,56.6,0);
  TGeoCombiTrans*   pMatrix26 = new TGeoCombiTrans("", dx,dy,dz,pMatrix27); 

  //Combi transformation:
  dx = 7.566000;
  dy = 0.000000;
  dz = sqrt(10.7*10.7-dx*dx);
  // Rotation:
  thx = -41.000000;    phx = 0.000000;
  thy = 90.000000;    phy = 90.000000;
  thz = 0.000000;    phz = 0.000000;
  //TGeoRotation *       pMatrix19 = new TGeoRotation("",thx,phx,thy,phy,thz,phz);
  TGeoRotation *       pMatrix19 = new TGeoRotation();
  pMatrix19->SetAngles(90,49,0);
  TGeoCombiTrans*   pMatrix18 = new TGeoCombiTrans("", dx,dy,dz,pMatrix19);

  //Combi transformation:
  dx = 7.566000;
  dy = 0.000000;
  dz = sqrt(10.7*10.7-dx*dx)+0.1;
  // Rotation:
  thx = -41.000000;    phx = 0.000000;
  thy = 90.000000;    phy = 90.000000;
  thz = 0.000000;    phz = 0.000000;
  //TGeoRotation *       pMatrix19 = new TGeoRotation("",thx,phx,thy,phy,thz,phz);
  TGeoRotation *       pMatrix23 = new TGeoRotation();
  pMatrix23->SetAngles(90,49,0);
  TGeoCombiTrans*   pMatrix22 = new TGeoCombiTrans("", dx,dy,dz,pMatrix23); 
	
  //Combi transformation:
  dx = 0.00000;
  dy = 0.00000;
  dz = 14.00000; //14
  // Rotation:
  thx = 90.000000;    phx = 0.000000;
  thy = 90.000000;    phy = 90.000000;
  thz = 0.000000;    phz = 0.000000;
  TGeoRotation *       pMatrix15 = new TGeoRotation("",thx,phx,thy,phy,thz,phz);
  TGeoCombiTrans*   pMatrix14 = new TGeoCombiTrans("", dx,dy,dz,pMatrix15); //Beam Tracker 1

  //Combi transformation:
  dx = 0.00000;
  dy = 0.00000;
  dz = 14.10000; //14
  // Rotation:
  thx = 90.000000;    phx = 0.000000;
  thy = 90.000000;    phy = 90.000000;
  thz = 0.000000;    phz = 0.000000;
  TGeoRotation *       pMatrix21 = new TGeoRotation("",thx,phx,thy,phy,thz,phz);
  TGeoCombiTrans*   pMatrix20 = new TGeoCombiTrans("", dx,dy,dz,pMatrix21); //Beam Tracker 2
	
  //Combi transformation:
  dx = 0.000000;
  dy = -2.100000;
  dz = 4.470000;
  // Rotation:
  thx = 0.000000;    phx = 0.000000;
  thy = 90.000000;    phy = 0.000000;
  thz = 90.000000;    phz = 90.000000;
  TGeoRotation *       pMatrix7 = new TGeoRotation("",thx,phx,thy,phy,thz,phz);
  TGeoCombiTrans*   pMatrix6 = new TGeoCombiTrans("", dx,dy,dz,pMatrix7); //Down
  
  //Combi transformation:
  dx = 0.000000;
  dy = 2.100000;
  dz = 4.470000;
  // Rotation:
  thx = 180.000000;    phx = 0.000000;
  thy = 90.000000;    phy = 0.000000;
  thz = 90.000000;    phz = 270.000000;
  TGeoRotation *	    pMatrix9 = new TGeoRotation("",thx,phx,thy,phy,thz,phz);
  TGeoCombiTrans*   pMatrix8 = new TGeoCombiTrans("", dx,dy,dz,pMatrix9); //Up
  
  //Combi transformation:
  dx = 2.100000;
  dy = 0.000000;
  dz = 4.470000;
  // Rotation:
  thx = 180.000000;    phx = 0.000000;
  thy = 90.000000;    phy = 90.000000;
  thz = 90.000000;    phz = 0.000000;
  TGeoRotation *        pMatrix11 = new TGeoRotation("",thx,phx,thy,phy,thz,phz);
  TGeoCombiTrans*   pMatrix10 = new TGeoCombiTrans("", dx,dy,dz,pMatrix11); //Left
  
  //Combi transformation:
  dx = -2.100000;
  dy = 0.000000;
  dz = 4.470000;
  // Rotation:
  thx = 0.000000;    phx = 0.000000;
  thy = 90.000000;    phy = 90.000000;
  thz = 90.000000;    phz = 180.000000;
  TGeoRotation *	    pMatrix13 = new TGeoRotation("",thx,phx,thy,phy,thz,phz);
  TGeoCombiTrans*   pMatrix12 = new TGeoCombiTrans("", dx,dy,dz,pMatrix13); //Right
 

  /************ Assembling everything together ****************/
  TGeoVolume *aTra = new TGeoVolumeAssembly("ATRA");
    
  //PCB
  aTra->AddNode(PCBvol1,11, pMatrix4);
  aTra->AddNode(PCBvol2,12, pMatrix2);
  aTra->AddNode(PCBvol2,13, pMatrix16);
  aTra->AddNode(PCBvol1,14, pMatrix18);
  aTra->AddNode(PCBvol1,15, pMatrix14); //Tracker
  aTra->AddNode(PCBvol3,16, pMatrix26);
  aTra->AddNode(PCBvol3,17, pMatrix28);
  aTra->AddNode(PCBvol3,18, pMatrix24);
  aTra->AddNode(PCBvol3,19, pMatrix22);
  aTra->AddNode(PCBvol3,10, pMatrix20); //Tracker

  //Si sensors
  aTra->AddNode(TraLog,1, pMatrix4);
  aTra->AddNode(TraLog,2, pMatrix2);
  aTra->AddNode(TraLog,3, pMatrix16);
  aTra->AddNode(TraLog,4, pMatrix18);
  aTra->AddNode(TraLog,5, pMatrix14); //Tracker
  aTra->AddNode(TraLog,6, pMatrix26);
  aTra->AddNode(TraLog,7, pMatrix28);
  aTra->AddNode(TraLog,8, pMatrix24);
  aTra->AddNode(TraLog,9, pMatrix22);
  aTra->AddNode(TraLog,10, pMatrix20); //Tracker
  
  // aTra->AddNode(TraLog,1, pMatrix6);//down
  // aTra->AddNode(TraLog,2, pMatrix8);//up
  // aTra->AddNode(TraLog,3, pMatrix10);//left
  // aTra->AddNode(TraLog,4, pMatrix12);//right
	
  
  
  TGeoRotation *rotg = new TGeoRotation();
  rotg->RotateX(0.);
  rotg->RotateY(0.);
  rotg->RotateZ(0.);
  dx=tx=0.0;
  dy=ty=0.0;
  dz=tz=0.0;
  
  TGeoCombiTrans *t_zero = new TGeoCombiTrans("t_zero");

  TGeoCombiTrans *t0 = new TGeoCombiTrans(tx,ty,tz,rotg);
  pWorld->AddNode(aTra, 1, GetGlobalPosition(t0));
   
  // ---------------   Finish   -----------------------------------------------
  gGeoMan->CloseGeometry();
  gGeoMan->CheckOverlaps(0.001);
  gGeoMan->PrintOverlaps();
  gGeoMan->Test();

  TFile* geoFile = new TFile(geoFileName, "RECREATE");
  top->Write();
  geoFile->Close();

  std::cout << "Creating geometry: "<<geoFileName<< std::endl;

  // --------------------------------------------------------------------------
}



TGeoCombiTrans* GetGlobalPosition(TGeoCombiTrans *fRef)
{
  if (fLocalTrans == kTRUE ) {
    
    if ( ( fThetaX == 0 )  && ( fThetaY==0 )  && ( fThetaZ == 0 )
	 &&
	 ( fX == 0 ) && ( fY == 0 ) && ( fZ == 0 )
	 )  return fRef;
    
    
    // X axis
    Double_t xAxis[3] = { 1. , 0. , 0. };
    Double_t yAxis[3] = { 0. , 1. , 0. };
    Double_t zAxis[3] = { 0. , 0. , 1. };
    // Reference Rotation
    fRefRot = fRef->GetRotation();
    
    if (fRefRot) {
      Double_t mX[3] = {0.,0.,0.};
      Double_t mY[3] = {0.,0.,0.};
      Double_t mZ[3] = {0.,0.,0.};
      
      fRefRot->LocalToMasterVect(xAxis,mX);
      fRefRot->LocalToMasterVect(yAxis,mY);
      fRefRot->LocalToMasterVect(zAxis,mZ);
      
      Double_t a[4]={ mX[0],mX[1],mX[2], fThetaX };
      Double_t b[4]={ mY[0],mY[1],mY[2], fThetaY };
      Double_t c[4]={ mZ[0],mZ[1],mZ[2], fThetaZ };
      
      ROOT::Math::AxisAngle aX(a,a+4);
      ROOT::Math::AxisAngle aY(b,b+4);
      ROOT::Math::AxisAngle aZ(c,c+4);
      
      ROOT::Math::Rotation3D fMatX( aX );
      ROOT::Math::Rotation3D fMatY( aY );
      ROOT::Math::Rotation3D fMatZ( aZ );
      
      ROOT::Math::Rotation3D  fRotXYZ = (fMatZ * (fMatY * fMatX));
      
      //cout << fRotXYZ << endl;
      
      Double_t fRotable[9]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
      fRotXYZ.GetComponents(
			    fRotable[0],fRotable[3],fRotable[6],
			    fRotable[1],fRotable[4],fRotable[7],
			    fRotable[2],fRotable[5],fRotable[8]
			    );
      TGeoRotation *pRot = new TGeoRotation();
      pRot->SetMatrix(fRotable);
      TGeoCombiTrans *pTmp = new TGeoCombiTrans(*fGlobalTrans,*pRot);
      
      // ne peut pas etre applique ici
      // il faut differencier trans et rot dans la multi.
      TGeoRotation rot_id;
      rot_id.SetAngles(0.0,0.0,0.0);
      
      TGeoCombiTrans c1;
      c1.SetRotation(rot_id);
      const Double_t *t = pTmp->GetTranslation();
      c1.SetTranslation(t[0],t[1],t[2]);
      
      TGeoCombiTrans c2;
      c2.SetRotation(rot_id);
      const Double_t *tt = fRefRot->GetTranslation();
      c2.SetTranslation(tt[0],tt[1],tt[2]);
      
      TGeoCombiTrans cc = c1 * c2 ;
      
      TGeoCombiTrans c3;
      c3.SetRotation(pTmp->GetRotation());
      TGeoCombiTrans c4;
      c4.SetRotation(fRefRot);
      
      TGeoCombiTrans ccc = c3 * c4;
      
      TGeoCombiTrans pGlobal;
      pGlobal.SetRotation(ccc.GetRotation());
      const Double_t *allt = cc.GetTranslation();
      pGlobal.SetTranslation(allt[0],allt[1],allt[2]);
      
      return  ( new TGeoCombiTrans( pGlobal ) );
      
    }else{
      
      cout << "-E- R3BDetector::GetGlobalPosition() \
	      No. Ref. Transformation defined ! " << endl;
      cout << "-E- R3BDetector::GetGlobalPosition() \
	      cannot create Local Transformation " << endl;
      return NULL;
    } //! fRefRot
    
  } else {
    // Lab Transf.
    if ( ( fPhi == 0 )  && ( fTheta==0 )  && ( fPsi == 0 )
	 &&
	 ( fX == 0 ) && ( fY == 0 ) && ( fZ == 0 )
	 )  return fRef;
    
    
    return ( new TGeoCombiTrans(*fGlobalTrans,*fGlobalRot) );
    
  }
}

