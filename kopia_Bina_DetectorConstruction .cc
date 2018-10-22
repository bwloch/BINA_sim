#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trap.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4Polyhedra.hh"
#include "G4Polycone.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4VSolid.hh"

#include "G4RunManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4UImanager.hh"

#include "Bina_DetectorConstruction.hh"
#include "Bina_DetectorMessenger.hh"
#include "Bina_Graniastoslupy.hh"

#include "Bina_PrimaryGeneratorAction.hh"

using namespace CLHEP;

Bina_DetectorConstruction::Bina_DetectorConstruction()
:
//npd_choice(0),

 Mars_sol(0),		Mars_log(0),		Mars_phs(0)
 ,Mars_vac_sol(0),	Mars_vac_log(0),	Mars_vac_phs(0)
 ,trgt_sol(0),		trgt_log(0),		trgt_phs(0)
 ,trgt_win_sol(0),	trgt_win_log(0)
 ,mwpc1_sol(0),		mwpc1_minus_sol(0),	mwpc1_log(0),		mwpc1_phs(0)
 ,mwpc1_frame_sol(0),	mwpc1_frame_log(0),	mwpc1_frame_phs(0)
 ,mwpc1_gas_sol(0),	mwpc1_gas_minus_sol(0),	mwpc1_gas_log(0),	mwpc1_gas_phs(0)
 ,mwpc1_wire_plane_mother_sol(0),	mwpc1_hole_sol(0),	mwpc1_wire_plane_mother_minus_sol(0)
 ,mwpc1_wire_sol(0),	mwpc1_wire_log(0)
 ,mwpc1_catode_foil_sol(0),mwpc1_catode_minus_sol(0),mwpc1_catode_foil_log(0)
 ,mwpc1_ring_sol(0),	mwpc1_ring_minus_sol(0),mwpc1_ring_log(0),	mwpc1_ring_phs(0)
 ,mwpc1_win_sol(0)
 ,pipe_sol(0),	pipe_log(0),	pipe_phs(0)
 ,win_scatering_sol(0),	win_scatering_log(0),	win_scatering_phs(0)
 //,fpMagField(0)
 ,walec(0)
 ,deltaE_w_sol(0)
//Bina E, dE
 ,Bina_sol(0),Bina_log(0),Bina_minus(0),Bina_phs(0)
 ,blok_sol(0)
//delta E
 ,deltae_W_sol(0)
 ,foil_sol(0), foil_min_sol(0), foil_log(0), foil_phs(0)
//Wall E
  ,wall_L_sol(0)  ,wall_sol(0)  ,wall_Spec_sol(0)  ,wall_C_sol(0)
 // ,wall_log()  ,wall_C_log(0)  ,wall_L_log(0)  //,wall_Spec_log(0)
 /* ,wall_phs() ,wallf_phs()*/


{
 // fpMagField = new Bina_MagneticField();
 int i;

  for (i=0;i<12;i++) { mwpc1_wire_plane_mother_log[i] = 0; mwpc1_wire_plane_mother_phs[i] = 0;}
  for (i=0;i<30;i++) { mwpc1_wire_c_sol[i] = 0;	mwpc1_wire_c_log[i] = 0; }
  for (i=0;i<150;i++){ mwpc1_wire_45_sol[i] = 0; mwpc1_wire_45_log[i] = 0; }
  for (i=0;i<250;i++){ mwpc1_wire_x_phs[i] = 0; mwpc1_wire_y_phs[i] = 0; }
  for (i=0;i<350;i++) mwpc1_wire_45_phs[i] = 0;
  for (i=0;i<8;i++)   mwpc1_catode_foil_phs[i] = 0;

//Bina E, dE  
  for (i=0;i<7;i++)
  {
   Edet_B_sol[i] = 0;    Edet_B_log[i] = 0;
   deltaE_B_sol[i] = 0;  deltaE_B_log[i] = 0;
  }
  for (i=0;i<200;i++)
  { Edet_B_phs[i]=0;     deltaE_B_phs[i]=0; } 
  
//delta E 
  for (i=0;i<10;i++)
  {
   deltae_Wc_sol[i]=0;   deltae_W_log[i]=0;
  } 
  for (i=0;i<40;i++)  { deltae_W_phs[i]=0;  }
//wall E
  for (i=0;i<10;i++)
  {
   wall_phs[i]=0;
  }
  
//default values for geometry
  #include "Bina_Detector.cfg"

  w_n_bar = &n_bar;
  w_n_pipe = &n_pipe;
  w_eps_separ = &eps_separ;
  w_z_shift = &z_shift;
  w_z_pipe = &z_pipe;


  detectorMessenger = new Bina_DetectorMessenger(this);
}


Bina_DetectorConstruction::~Bina_DetectorConstruction()
{
delete detectorMessenger;
}

void Bina_DetectorConstruction::Seen()
{
  Mars_log->SetVisAttributes (G4VisAttributes::Invisible);
  mwpc1_log->SetVisAttributes (G4VisAttributes::Invisible);
  for (int i=0;i<9;i++) mwpc1_wire_plane_mother_log[i]->SetVisAttributes (G4VisAttributes::Invisible);
}


G4VPhysicalVolume* Bina_DetectorConstruction::Construct()
{

  #include "Bina_Material.cfg"
//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------

  double inRadius, outRadius, hight, startAngle, endAngle;
  G4double box_x, box_y, box_z, box_z_vac;
  double temp,temp2;
  G4ThreeVector position;
  const double pi = 3.1415926535;
  int m,i;
  G4RotationMatrix* rot1[5];
  G4RotationMatrix* rot2[5];
  double ang_dif,ang;


  G4double pipe_dim[21] = {  0.0,  3.56, 3.75,
     		 	3.5,  3.60, 3.80,
     		 	10.2, 3.60, 3.80,	//!!!!!!!! uwaga zmiana ostatniej wartosci do 3.80 z 4.60
     		 	13.0, 4.30, 4.60,
			39.7, 4.30, 4.60,	//zmiana z 4.80 do 4.60 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			40.6, 4.50, 4.80,
			63.6, 4.50, 4.80};
//  G4double pipe_lock[2] = {1.,2.};

  //------------------------------
  //		World
  //------------------------------

  box_x = dim_mars[0];
  box_y = dim_mars[1];
  box_z = dim_mars[2];

  Mars_sol= new G4Box("Mars_world_solid",box_x*cm,box_y*cm,box_z*cm);
  Mars_log= new G4LogicalVolume(Mars_sol, Air, "Mars_world_logical");

  box_x = dim_mars[0]-10*eps_separ;
  box_y = dim_mars[1]-10*eps_separ;
  box_z_vac = (dim_mars[2] + exit_win_z - exit_win_th)/2.-10*eps_separ ;
  z_shift = box_z - box_z_vac - 10*eps_separ; // shift of positions for objects
                                              // placed in vacuum

  Mars_vac_sol= new G4Box("Mars_vacuum_solid",box_x*cm,box_y*cm,box_z_vac*cm);
  Mars_vac_log= new G4LogicalVolume(Mars_vac_sol,Vacuum,"Mars_vacuum_logical");

  //------------------------------
  //	       Beam pipe 
  //------------------------------

  // - > steal tube

  double phiStart = 0.;			// phi_min
  double phiTotal = 360.;		// phi_max
  int nZplanes = n_pipe + 1;		//  number of segments + 1

  G4double* Zplanes = new G4double[n_pipe + 1];	//z coordinate - starting position
  G4double* r_in = new G4double[n_pipe + 1];	//r_in coordinate of these corners
  G4double* r_out = new G4double[n_pipe + 1];	//r_out coordinate of these corners

  Zplanes[0] = pipe_dim[0]*cm;
  r_in[0] = 0.*cm;
  r_out[0] = pipe_dim[2]*cm;

  for (i=1;i<n_pipe+1;i++)
  {
    Zplanes[i] = (Zplanes[0] + pipe_dim[(i-1)*3 + 3])*cm;
    r_in[i] = 0.*cm;
    r_out[i] = (pipe_dim[(i-1)*3 + 5])*cm;
  }

  pipe_sol = new G4Polycone("Pipe_solid",phiStart*deg,phiTotal*deg,nZplanes,Zplanes,r_in,r_out);
  pipe_log = new G4LogicalVolume(pipe_sol, SSteel, "Pipe_logical");

  // - > vacum

  r_out[0] = pipe_dim[1]*cm;

  // -> lock

  //inRadius = pipe_dim[2];
  //outRadius = inRadius + pipe_lock[0];
  //hight = pipe_lock[1];
  //startAngle = 0.;
  //endAngle = 360.;

  //------------------------------
  // MWPC1 (CF4+Isobutane) with vacuum hole
  //------------------------------

  // - > Mother volume (air)
  box_x = 1.1*mwpcdimfx;
  box_y = 1.1*mwpcdimfx;
  box_z = 1.1*mwpcdimfz;

  position = G4ThreeVector(0.,0.,- mwpcZplace);
  mwpc1_sol = new G4Box("Mwpc1_solid_mother",box_x,box_y,box_z);
  mwpc1_minus_sol = new G4SubtractionSolid("Mwpc1_logical_minus_pipe_mother",mwpc1_sol,pipe_sol,0,position);
  mwpc1_log = new G4LogicalVolume(mwpc1_minus_sol, Air, "Mwpc1_logical_mother", 0, 0, 0);
  // - > frame (Al)

  phiStart = 45.;
  phiTotal = 360.;
  int nSide = 4;		//number of segments
  nZplanes = 2;			//number corners    ??? rs

  const G4double Zplanes_1[2] = {-mwpcdimfz, mwpcdimfz};  //z coordinate of these corners
  const G4double r_in_1[2] =    {mwpcdimgx, mwpcdimgx };  //r_in these corners
  const G4double r_out_1[2] =   {mwpcdimfx, mwpcdimfx};  //r_out

  mwpc1_frame_sol = new G4Polyhedra("Mwpc1_frame_solid",phiStart*deg,phiTotal*deg,nSide,nZplanes,Zplanes_1,r_in_1,r_out_1);
  mwpc1_frame_log = new G4LogicalVolume(mwpc1_frame_sol, Al, "Mwpc1_frame_logical");

  // - > gas_chamber

  box_x = mwpcdimgx;
  box_y = mwpcdimgx;
  box_z = mwpcdimgz;

  position = G4ThreeVector(0.,0.,- mwpcZplace);
  mwpc1_gas_sol = new G4Box("Mwpc1_gas_solid",box_x,box_y,box_z);
  mwpc1_gas_minus_sol = new G4SubtractionSolid("Mwpc1_gas_minus_pipe_logical",mwpc1_gas_sol,pipe_sol,0,position);
  mwpc1_gas_log = new G4LogicalVolume(mwpc1_gas_minus_sol, GasMWPC, "Mwpc1_gas_logical");

  // -> ring

  inRadius = mwpcholeIn; //dim_ring[0];
  outRadius = mwpcholeOut; //dim_ring[1];
  hight = mwpcdimgz; //dim_mwpc1[3];
  startAngle = 0.;
  endAngle = 360.;

  position = G4ThreeVector(0.,0.,- mwpcZplace);
  mwpc1_ring_sol = new G4Tubs("Mwpc1_ring_solid",inRadius,outRadius,hight,startAngle*deg,endAngle*deg);
  mwpc1_ring_minus_sol = new G4SubtractionSolid("Mwpc1_ring_minus_pipe_logical",mwpc1_ring_sol,pipe_sol,0,position);
  mwpc1_ring_log = new G4LogicalVolume(mwpc1_ring_minus_sol, Al, "Mwpc1_ring_logical");

  // - > Mwpc1 window (Kapton)

  box_x = mwpcdimgx; 
  box_y = mwpcdimgx; 
  box_z = par_mwpc[0]*cm;

  temp = -mwpcdimgz - par_mwpc[0]*cm;
  position = G4ThreeVector(0.,0.,z_pipe*cm - temp);
  mwpc1_win_sol = new G4Box("Mwpc1_window_solid",box_x,box_y,box_z);
  mwpc1_win_minus_sol[0] = new G4SubtractionSolid("Mwpc1_window_minus_pipe_left_logical",mwpc1_win_sol,pipe_sol,0,position);
  mwpc1_win_log[0] = new G4LogicalVolume(mwpc1_win_minus_sol[0], Kapton, "Mwpc1_window_left_logical");

  position = G4ThreeVector(0.*cm,0.*cm, z_pipe*cm + temp);
  mwpc1_win_minus_sol[1] = new G4SubtractionSolid("Mwpc1_window_minus_pipe_right_logical",mwpc1_win_sol,pipe_sol,0,position);
  mwpc1_win_log[1] = new G4LogicalVolume(mwpc1_win_minus_sol[1], Kapton, "Mwpc1_window_right_logical");

  //divide chamber wire

  double r_wire = par_mwpc[1];     // 10 um
  double wire_space = par_mwpc[4];   // 2 mm
  double catode_space = par_mwpc[3]; // 4 mm *0.5 (half dimension)
  double start_pos = par_mwpc[4];
  double start_wire = par_mwpc[2];


  // mother wire plane

  box_x =  mwpcdimgx;
  box_y =  mwpcdimgx;
  box_z = 0.0013*cm;//catode_space;	// r_wire < catode_foil < wire_mother
 // box_z = 0.0003*cm;//catode_space;	// r_wire < catode_foil < wire_mother
  mwpc1_wire_plane_mother_sol = new G4Box("Mwpc1_wire_plane_mother_solid",box_x,box_y,box_z);

  inRadius = 0.;
  outRadius =  mwpcholeOut;
  hight = 3.*catode_space*cm;
  startAngle = 0.;
  endAngle = 360.;

  mwpc1_hole_sol = new G4Tubs("Mwpc1_hole_solid",inRadius,outRadius,hight,startAngle*deg,endAngle*deg);
  mwpc1_wire_plane_mother_minus_sol = new G4SubtractionSolid("Mwpc1_mother_plane_minus_hole_logical"
  				,mwpc1_win_sol,mwpc1_hole_sol,0,position);
  for (m=0;m<10;m++) mwpc1_wire_plane_mother_log[m] = new G4LogicalVolume(mwpc1_wire_plane_mother_sol, GasMWPC,"Mwpc1_wire_plane_mother_logical");

  box_z = 0.00125*cm;
  mwpc1_catode_foil_sol = new G4Box("Mwpc1_catode_solid",box_x,box_y,box_z);
  mwpc1_catode_foil_log = new G4LogicalVolume(mwpc1_catode_foil_sol,Mylar,"Mwpc1_catode_foil_logical");

 
  walec = new G4Tubs("walec",0.*cm,5.*cm,20.*cm,0.,2*M_PI);



//--------- Position Physical Volumes ---------

  //------------------------------
  //		World
  //------------------------------

  //  Must place the World Physical volume unrotated at (0,0,0).

  Mars_phs = new G4PVPlacement(0,               // no rotation
                                 G4ThreeVector(), // at (0,0,0)
				 "Mars_world_physical",         // its name
                                 Mars_log,      // its logical volume
                                 0,               // its mother  volume (Physical)
                                 false,           // no boolean operations
                                 0);              // no field specific to volume
  Mars_vac_phs = new G4PVPlacement(0,               // no rotation
                                 G4ThreeVector(0,0,-z_shift*cm), 
				 "Mars_vacuum_physical",         // its name
                                 Mars_vac_log,      // its logical volume
                                 Mars_phs,               // its mother  volume (Physical)
                                 false,           // no boolean operations
                                 0);              // no field specific to volume

 
  //------------------------------
  //		Target
  //------------------------------

  if(targetis)
  {

    // Target frame
     //  trgt_fr_sol = new G4Tubs("Tframe_solid",toutrad+eps_separ*cm,tfroutrad,thigh,0.,2.*M_PI);
  trgt_fr_sol = new G4Tubs("Tframe_solid",0,tfroutrad,thigh,0.,2.*M_PI);
  trgt_fr_log = new G4LogicalVolume(trgt_fr_sol, SSteel, "Tframe_logical");
  
  // - > liquid hydrogen
  trgt_sol = new G4Tubs("Target_liquid_solid",0.,toutrad,thigh,0.,2.*M_PI);
  trgt_log = new G4LogicalVolume(trgt_sol, lD2, "Target_liquid_logical");

  // - > target window (aramica 25 um); here in [mm]

  hight = 0.000000001;

  trgt_win_sol = new G4Tubs("Target_window_solid",0.,toutrad,hight,0.,2.*M_PI);
  trgt_win_log = new G4LogicalVolume(trgt_win_sol, Aramica, "Target_window_logical");

  // - > target Al shielding
  trgt_sh_sol = new G4Tubs("Tframe_solid",4.95*cm,5.*cm,10*cm,0.,2.*M_PI);
  walec = new G4Tubs("walec",0.,8*cm,5.*cm,0.,2*M_PI);
  G4RotationMatrix* rothole = new G4RotationMatrix();
  rothole->rotateX(M_PI/2.);
  position = G4ThreeVector(0,-5*cm,0);
  trgt_shc_sol = new G4SubtractionSolid("tr_shc_sol",
       trgt_sh_sol,walec,rothole,position);
  trgt_sh_log = new G4LogicalVolume(trgt_shc_sol, SSteel, "Tshield_logical");



  //------------------------------
  //		Target
  //------------------------------
  
  // - > liquid hydrogen

  position = G4ThreeVector(tXplace, tYplace,tZplace + z_shift*cm);
  G4RotationMatrix* rotshield = new G4RotationMatrix();
  rotshield->rotateX(-M_PI/2.);
  trgt_fr_phs = new G4PVPlacement(0,position,"Tframe_physical",trgt_fr_log,Mars_vac_phs,false,0);
  trgt_phs = new G4PVPlacement(0,G4ThreeVector(0,0,0),"Target_physical",trgt_log,trgt_fr_phs,false,0);
  //  trgt_sh_phs = new G4PVPlacement(rotshield,position,"Tshield_physical",trgt_sh_log,Mars_vac_phs,false,0);

  // - > target window (aramica 25 um)
  temp = tZplace - ( 0.00125 + eps_separ)*cm -thigh;
  position = G4ThreeVector(tXplace, tYplace, temp+z_shift*cm);
  trgt_win_phs[0] = new G4PVPlacement(0,position,"Target_window1_physical",trgt_win_log, Mars_vac_phs,false,1);

  temp = tZplace + (0.00125 + eps_separ)*cm + thigh;
  position = G4ThreeVector(tXplace, tYplace ,temp+z_shift*cm);
  trgt_win_phs[1] = new G4PVPlacement(0,position,"Target_window2_physical",trgt_win_log,Mars_vac_phs,false,2);
  }

  ///////////////////////
  ///  Exit window of scattering chamber
  ///////////////////////
  
  exit_win_sol=new G4Tubs("exit_win_sol",0.,exit_win_rad*cm,exit_win_th*cm,0.,2*M_PI);
  position = G4ThreeVector(0.,0.,(exit_win_z+exit_win_th)*cm); 
  exit_win_log = new G4LogicalVolume(exit_win_sol, Aramica, "exit_win_log"); 
  exit_win_phs = new G4PVPlacement(0,position,"Exit_win_phs",exit_win_log,Mars_phs,false,0);

  //    G4VisAttributes* kolorf= new G4VisAttributes(G4Colour(0.,1.0,1.));
  //  exit_win_log ->SetVisAttributes(kolorf);


  //------------------------------
  // MWPC1 (CF4+Isobutane) with vacuum hole
  //------------------------------
 if(mwpcis)
 {
  {
  // - > Mother volume (air)

  position = G4ThreeVector(mwpcXplace,mwpcYplace,mwpcZplace);
  mwpc1_phs = new G4PVPlacement(0,position,"Mwpc1_physical",mwpc1_log,Mars_phs,false,0);

  // - > frame (Al)
  position = G4ThreeVector();
  mwpc1_frame_phs = new G4PVPlacement(0,position,"Mwpc1_frame_physical",mwpc1_frame_log,mwpc1_phs,false,0);

  // - > gas
  position = G4ThreeVector();
  mwpc1_gas_phs = new G4PVPlacement(0,position,"Mwpc1_gas_physical",mwpc1_gas_log,mwpc1_phs,false,0);

  // - > Mwpc1 windows ???????
  temp = -mwpcdimgz - par_mwpc[0]*cm; //change from 1 (3)
  position = G4ThreeVector(0., 0., temp);
  mwpc1_win_phs[0] = new G4PVPlacement(0,position,"Mwpc1_window_left_physical",mwpc1_win_log[0],mwpc1_phs,false,1);
  mwpc1_win_phs[1] = new G4PVPlacement(0,-position,"Mwpc1_window_right_physical",mwpc1_win_log[1],mwpc1_phs,false,2);

  // - > vacuum hole

  position = G4ThreeVector();
  mwpc1_ring_phs = new G4PVPlacement(0,position,"Mwpc1_ring_center_physical",mwpc1_ring_log,mwpc1_gas_phs,false,0);

  // - > wire plane + catode foil

  inRadius = 0.;
  outRadius = r_wire*cm;
  hight = mwpcdimgx;
  startAngle = 0.;
  endAngle = 360.;

  mwpc1_wire_sol = new G4Tubs("Mwpc1_wire_solid",inRadius,outRadius,hight,startAngle*deg,endAngle*deg);
  mwpc1_wire_log = new G4LogicalVolume(mwpc1_wire_sol, W, "Mwpc1_wire_logical");
  mwpc1_wire_log->SetVisAttributes (G4VisAttributes::Invisible);

  double max_size =  mwpcdimgz - catode_space*cm - start_pos*cm;
  double max_size2 = mwpcdimgx - start_wire*cm;
  double max_size3 = sqrt(2.)* mwpcdimgx - 6.*wire_space*cm - start_wire*cm;
  temp = - mwpcdimgz + start_pos*cm + catode_space*cm;
  temp2 = - mwpcdimgx + start_wire*cm + wire_space*cm;
  double delta;
  G4RotationMatrix* plane_x = new G4RotationMatrix(0.*deg,90.*deg,0.*deg);
  G4RotationMatrix* plane_y = new G4RotationMatrix(90.*deg,90.*deg,0.*deg);
  G4RotationMatrix* plane_45 = new G4RotationMatrix(45.*deg,90.*deg,0.*deg);

  i = -1;
  m = 0;
  int j=0,k=0,l=0,n=0,plane=1;

  temp2 = 0.5*wire_space*cm;
  while ((delta = sqrt(0.25*mwpcholeOut*mwpcholeOut - 0.25*temp2*temp2)) >= 0.)
  {
    hight = 0.5*mwpcdimgx - delta;
    mwpc1_wire_c_sol[k] = new G4Tubs("Mwpc1_wire_center_solid",inRadius,outRadius,hight,startAngle*deg,endAngle*deg);
    mwpc1_wire_c_log[k] = new G4LogicalVolume(mwpc1_wire_c_sol[k], W, "Mwpc1_wire_center_logical");
    mwpc1_wire_c_log[k]->SetVisAttributes (G4VisAttributes::Invisible);
    k++;
    temp2 += wire_space*cm;
  }

  k = 0;
  temp2 = 0.5*wire_space*cm;
  while ((delta = sqrt(0.25*mwpcholeOut*mwpcholeOut - 0.25*temp2*temp2)) >= 0.) //center part plane 45
  {
    hight = sqrt(2.)*0.5*mwpcdimgx - delta - 0.5*temp2;
    mwpc1_wire_45_sol[k] = new G4Tubs("Mwpc1_wire_45_center_solid",inRadius,outRadius,hight,startAngle*deg,endAngle*deg);
    mwpc1_wire_45_log[k] = new G4LogicalVolume(mwpc1_wire_45_sol[k], W, "Mwpc1_wire_45_center_logical");
    mwpc1_wire_45_log[k]->SetVisAttributes (G4VisAttributes::Invisible);
    k++;
    temp2 += wire_space*cm;
  }

  while(temp2 <= max_size3)	// rest part palne 45
  {
    hight = sqrt(2.)*mwpcdimgx - temp2;
    mwpc1_wire_45_sol[k] = new G4Tubs("Mwpc1_wire_45_solid",inRadius,outRadius,hight,startAngle*deg,endAngle*deg);
    mwpc1_wire_45_log[k] = new G4LogicalVolume(mwpc1_wire_45_sol[k], W, "Mwpc1_wire_45_logical");
    mwpc1_wire_45_log[k]->SetVisAttributes (G4VisAttributes::Invisible);
    temp2 += wire_space*cm;
    k++;
  }

  j = 0;

  do
  {
    i++;
    position = G4ThreeVector(0.,0.,temp);
    mwpc1_wire_plane_mother_phs[i] = new G4PVPlacement(0,position,"Mwpc1_slice_physical"
 						    ,mwpc1_wire_plane_mother_log[i],mwpc1_gas_phs,false,i+1);
    if (i%3 == 1) // anode
    {
      temp2 = 0.5*wire_space*cm;
      k = 0;
      if (plane == 1)	//vertical plane
      {
        while ((delta = sqrt(0.25*mwpcholeOut*mwpcholeOut - 0.25*temp2*temp2)) >= 0.)
        {
          position = G4ThreeVector(temp2, 0.5*mwpcdimgx + delta,0.);
          mwpc1_wire_x_phs[j] = new G4PVPlacement(plane_x,position,"Mwpc1_wire_x_c_u_physical",mwpc1_wire_c_log[k],mwpc1_wire_plane_mother_phs[i],false,j); ++j;
          mwpc1_wire_x_phs[j] = new G4PVPlacement(plane_x,-position,"Mwpc1_wire_x_c_u_physical",mwpc1_wire_c_log[k],mwpc1_wire_plane_mother_phs[i],false,j); ++j;
          position = G4ThreeVector(temp2,-(0.5*mwpcdimgx + delta),0.);
	  mwpc1_wire_x_phs[j] = new G4PVPlacement(plane_x,position,"Mwpc1_wire_x_c_u_physical",mwpc1_wire_c_log[k],mwpc1_wire_plane_mother_phs[i],false,j); ++j;
          mwpc1_wire_x_phs[j] = new G4PVPlacement(plane_x,-position,"Mwpc1_wire_x_c_u_physical",mwpc1_wire_c_log[k],mwpc1_wire_plane_mother_phs[i],false,j); ++j;
	  temp2 += wire_space*cm;
	  k++;
        }

        do
        {
          position = G4ThreeVector(temp2,0.,0.);
          mwpc1_wire_x_phs[j] = new G4PVPlacement(plane_x,position,"Mwpc1_wire_x_physical"
  	  					    ,mwpc1_wire_log,mwpc1_wire_plane_mother_phs[i],false,j); ++j;
          mwpc1_wire_x_phs[j] = new G4PVPlacement(plane_x,-position,"Mwpc1_wire_x_physical"
  						    ,mwpc1_wire_log,mwpc1_wire_plane_mother_phs[i],false,j); ++j;
	  temp2 += wire_space*cm;
        }
        while (temp2 <= max_size2);
      }

      else if (plane == 2)		//horizontal plane
      {
        while ((delta = sqrt(0.25*mwpcholeOut*mwpcholeOut - 0.25*temp2*temp2)) >= 0.)
        {
	  position = G4ThreeVector( 0.5*mwpcdimgx + delta,temp2,0.);
          mwpc1_wire_y_phs[m] = new G4PVPlacement(plane_y,position,"Mwpc1_wire_y_c_u_physical"
   					    ,mwpc1_wire_c_log[k],mwpc1_wire_plane_mother_phs[i],false,m); ++m;

          mwpc1_wire_y_phs[m] = new G4PVPlacement(plane_y,-position,"Mwpc1_wire_y_c_u_physical"
   					    ,mwpc1_wire_c_log[k],mwpc1_wire_plane_mother_phs[i],false,m); ++m;
          position = G4ThreeVector(-(0.5*mwpcdimgx + delta),temp2,0.);
	  mwpc1_wire_y_phs[m] = new G4PVPlacement(plane_y,position,"Mwpc1_wire_y_c_u_physical"
  					    ,mwpc1_wire_c_log[k],mwpc1_wire_plane_mother_phs[i],false,m); ++m;
          mwpc1_wire_y_phs[m] = new G4PVPlacement(plane_y,-position,"Mwpc1_wire_y_c_u_physical"
   					    ,mwpc1_wire_c_log[k],mwpc1_wire_plane_mother_phs[i],false,m); ++m;
	  temp2 += wire_space*cm;
	  k++;
        }

        do
        {
          position = G4ThreeVector(0.,temp2,0.);
          mwpc1_wire_y_phs[m] = new G4PVPlacement(plane_y,position,"Mwpc1_wire_y_physical"
  						    ,mwpc1_wire_log,mwpc1_wire_plane_mother_phs[i],false,m); ++m;
          mwpc1_wire_y_phs[m] = new G4PVPlacement(plane_y,-position,"Mwpc1_wire_y_physical"
  						    ,mwpc1_wire_log,mwpc1_wire_plane_mother_phs[i],false,m); ++m;
	  temp2 += wire_space*cm;
        }
        while (temp2 <= max_size2);
      }

      else // plane 45 deg
      {
        while ((delta = sqrt(0.25*mwpcholeOut*mwpcholeOut - 0.25*temp2*temp2)) >= 0.)
        {
	  hight = sqrt(2.)*0.5*mwpcdimgx - delta - 0.5*temp2;
          position = G4ThreeVector((cos(pi/4.)*(-hight - 2.*delta + temp2)),(sin(pi/4.)*(hight + 2.*delta + temp2)),0.);
          mwpc1_wire_45_phs[n] = new G4PVPlacement(plane_45,position,"Mwpc1_wire_45_physical"
   					    ,mwpc1_wire_45_log[k],mwpc1_wire_plane_mother_phs[i],false,n); ++n;

          mwpc1_wire_45_phs[n] = new G4PVPlacement(plane_45,-position,"Mwpc1_wire_45_physical"
   					    ,mwpc1_wire_45_log[k],mwpc1_wire_plane_mother_phs[i],false,n); ++n;
          position = G4ThreeVector((cos(pi/4.)*(hight + 2.*delta + temp2)),(sin(pi/4.)*(-hight - 2.*delta + temp2)),0.);
	  mwpc1_wire_45_phs[n] = new G4PVPlacement(plane_45,position,"Mwpc1_wire_45_physical"
   					    ,mwpc1_wire_45_log[k],mwpc1_wire_plane_mother_phs[i],false,n); ++n;
          mwpc1_wire_45_phs[n] = new G4PVPlacement(plane_45,-position,"Mwpc1_wire_45_physical"
   					    ,mwpc1_wire_45_log[k],mwpc1_wire_plane_mother_phs[i],false,n); ++n;
	  temp2 += wire_space*cm;
	  k++;
        }

        do
        {
          position = G4ThreeVector(cos(pi/4.)*temp2,sin(pi/4.)*temp2,0.);
          mwpc1_wire_45_phs[n] = new G4PVPlacement(plane_45,position,"Mwpc1_wire_x_physical"
  						    ,mwpc1_wire_45_log[k],mwpc1_wire_plane_mother_phs[i],false,n); ++n;
          mwpc1_wire_45_phs[n] = new G4PVPlacement(plane_45,-position,"Mwpc1_wire_x_physical"
  						    ,mwpc1_wire_45_log[k],mwpc1_wire_plane_mother_phs[i],false,n); ++n;
	  temp2 += wire_space*cm;
	  k++;
        }
        while (temp2 <= max_size3);
      };

      if (plane == 3) plane = 0;
      plane ++;
    }
    //catode
    else mwpc1_catode_foil_phs[l] =  new G4PVPlacement(0,G4ThreeVector(),"Mwpc1_catode_physical"
  						    ,mwpc1_catode_foil_log,mwpc1_wire_plane_mother_phs[i],false,l); ++l;

    temp += 2.*catode_space*cm;
  }
  while(temp <= max_size);
 } 
}
//////////////////////
/// Bina
//////////////////////
//----------------------------
//              Kulka
//----------------------------
 if(ballis>0)
 {
  
   G4double a1_5=4.357;
   G4double h1_5=3.067;
   G4double a2_5=3.735;
   G4double ak_5=4.3;
   G4double h_5=1.5;

   wymiary G5krotkie(a1_5,h1_5,a2_5,h_5,ak_5);
   G4double a2_52=2.488;
   G4double h_52=4.5;

   wymiary G5dlogie(a1_5,h1_5,a2_52,h_52,ak_5);
   G4double a1_6=4.357;
   G4double h1_6=3.885;
   G4double a2_6=3.735;
   G4double ak_6=7.2;
   G4double h_6=1.5;

   wymiary G6krotkie(a1_6,h1_6,a2_6,h_6,ak_6);
   G4double  a2_62=2.488;
   G4double h_62=4.5;

   wymiary G6dlogie(a1_6,h1_6,a2_62,h_62,ak_6);
   
   G4double R[8];
   R[0]=G5krotkie.get_r();
   R[1]=G5krotkie.get_rm();
   R[2]=G5dlogie.get_r();
   R[3]=G5dlogie.get_rm();
   R[4]=G6krotkie.get_r();
   R[5]=G6krotkie.get_rm();
   R[6]=G6dlogie.get_r();
   R[7]=G6dlogie.get_rm();
  
   Edet_B_sol[0]=przecinka(blok_sol,G5krotkie,1);
   Edet_B_log[0]= new G4LogicalVolume(Edet_B_sol[0],BC408,"E_short5_log");// E krotkieG5;
   deltaE_B_sol[0]=przecinka(blok_sol,G5krotkie,2);
   deltaE_B_log[0]= new G4LogicalVolume(deltaE_B_sol[0],BC444,"dE_short5_log");// dE krotkie G5;
   Edet_B_sol[1]=przecinka(blok_sol,G5dlogie,1);
   Edet_B_log[1]= new G4LogicalVolume(Edet_B_sol[1],BC408,"E_long5_log");// E dlogie  G5;
   deltaE_B_sol[1]=przecinka(blok_sol,G5dlogie,2);
   deltaE_B_log[1]= new G4LogicalVolume(deltaE_B_sol[1],BC444,"dE_long5_log");//  dloge G521 dE;
   Edet_B_sol[2]=przecinka(blok_sol,G6krotkie,1);
   Edet_B_log[2]= new G4LogicalVolume(Edet_B_sol[2],BC408,"E_short6_log");// E G6 krotki;
   deltaE_B_sol[2]=przecinka(blok_sol,G6krotkie,2);
   deltaE_B_log[2]= new G4LogicalVolume(deltaE_B_sol[2],BC444,"dE_short6_log");// krotki  dE G611;
   Edet_B_sol[3]=przecinka(blok_sol,G6dlogie,1);
   Edet_B_log[3]= new G4LogicalVolume(Edet_B_sol[3],BC408,"E_long6_log");// dlogie G62 E;
   deltaE_B_sol[3]=przecinka(blok_sol,G6dlogie,2);
   deltaE_B_log[3]= new G4LogicalVolume(deltaE_B_sol[3],BC444,"dE_long6_log");// dlogie G621 dE;
   Edet_B_sol[4]=przecinka(blok_sol,G6dlogie,1,1);
   Edet_B_log[4]= new G4LogicalVolume(Edet_B_sol[4],BC408,"E_long6_spec1_log");// G1 E;
   deltaE_B_sol[4]=przecinka(blok_sol,G6dlogie,2,1);
   deltaE_B_log[4]= new G4LogicalVolume(deltaE_B_sol[4],BC444,"dE_long6_spec1_log");// G11 dE;
   Edet_B_sol[5]=przecinka(blok_sol,G6dlogie,1,2);
   Edet_B_log[5]= new G4LogicalVolume(Edet_B_sol[5],BC408,"E_long6_spec2_log");// G2 E;
   deltaE_B_sol[5]=przecinka(blok_sol,G6dlogie,2,2);
   deltaE_B_log[5]= new G4LogicalVolume(deltaE_B_sol[5],BC444,"dE_long6_spec2_log");// G21 dE;
   Edet_B_sol[6]=przecinka(blok_sol,G6dlogie,1,3);
   Edet_B_log[6]= new G4LogicalVolume(Edet_B_sol[6],BC408,"E_long6_spec3_log");// G3 E;
   deltaE_B_sol[6]=przecinka(blok_sol,G6dlogie,2,3);
   deltaE_B_log[6]= new G4LogicalVolume(deltaE_B_sol[6],BC444,"dE_long6_spec3_log");// G31 dE;
   

   int NumBina[149]={136,137,113,112,111,135, 
		     148,149,133,132,131,147, 
		     145,146,128,127,126,144, 
		     142,143,123,122,121,141, 
		     139,140,118,117,116,138,

		     108,107,134,110,109,  
		     102,101,129,130,103,
		     96,95,124,125,97,  
		     90,89,119,120, 91,
		     84,83,114,115,85,
		  
		     105,106,79,78,77,104,
		     99,100,73,72,71,98,   
		     93,94,67,66,65,92,  
		     87,88,61,60,59,86,
		     
		     81,82,52,51,50,80,
		     75,76,46,45,44,74,  
		     69,70,40,39,38,68, 
		     63,64,34,33,32,62,   
		     57,58,28,27,26,56,
		  
		     54,55,25,24,53,  
		     48,49,20,19,47,   
		     42,43,15,14,41,  
		     36,37,10,9,35,     
		     30,31,5,4,29,
		  
		     22,23,21, 
		     17,18,16,  
		     12,13,11,     
		     7,8,6,  
		     2,3,1};

   int n_el=350;
   polozenia* plac=new polozenia[400];
   planowanie(plac,R,n_el);

    /////////////////////////
   int gE = 0, gdE=0, nrde=0 ; 
   G4double xb, yb, zb;
   G4RotationMatrix* rotBin = new G4RotationMatrix[400];
   for(int t=0;t<n_el;t++)//n_el  242-270 dobrze widac, jakby obrocone bylo
     {
       if(plac[t].typ!=-1)
	 {
//diagnostyka:::	 G4cout<<t<</*' '<<plac[t].rx<<' '<<plac[t].ry<<' '<<plac[t].rz<<' '
//	 <<plac[t].x<<' '<<plac[t].y<<' '<<plac[t].z<<' '<<plac[t].typ<<
//	 "\n    "<<xb<<' '<<yb<<' '<<zb<<' '<<plac[t].rz<<"=phi, cosphi="<<cos(plac[t].rz)<<' '<<plac[t].ry<<"=phi, cosphi="<<cos(plac[t].ry)<<'\n';
//do{	   if(plac[t].rx>2*M_PI){plac[t].rx=plac[t].rx-2*M_PI;}
//	   if(plac[t].ry>2*M_PI){plac[t].ry=plac[t].ry-2*M_PI;}
//	   if(plac[t].z>2*M_PI){plac[t].z=plac[t].z-2*M_PI;}}
//	   while (plac[t].rx>2*M_PI||plac[t].ry>2*M_PI||plac[t].z>2*M_PI);
	   //??????
	   xb =  plac[t].x*cm *sin(plac[t].y)*cos(plac[t].z) + ballXplace;
	   yb =  plac[t].x*cm *sin(plac[t].y)*sin(plac[t].z) + ballYplace;
	   zb =  plac[t].x*cm *cos(plac[t].y) + ballZplace + z_shift*cm;
	   position.setRThetaPhi ((plac[t].x*cm),plac[t].y,plac[t].z);
	   rotBin->rotateX(plac[t].rx);//1.10724);
	   rotBin->setPhi(plac[t].rz);
	   rotBin->setPsi(plac[t].ry); //!
	   if(plac[t].typ==0||plac[t].typ==1||plac[t].typ==2||plac[t].typ==3|| plac[t].typ==4||plac[t].typ==5||plac[t].typ==6)
	     { 
              //Edet_B_phs[gE]=new G4PVPlacement(rotBin,position,"EDet_BINA", Edet_B_log[plac[t].typ],  Mars_phs,false,NumBina[gE]); gE++;
              Edet_B_phs[gE]=new G4PVPlacement(rotBin, G4ThreeVector(xb, yb, zb),"EDet_BINA",  Edet_B_log[plac[t].typ],  Mars_vac_phs,false,NumBina[gE]); gE++;
	     
	       }
	   if(plac[t].typ==10||plac[t].typ==11||plac[t].typ==12||plac[t].typ==13|| plac[t].typ==14||plac[t].typ==15||plac[t].typ==16)
	     { 
	       nrde = plac[t].typ-10;
	       //deltaE_B_phs[gdE]=new G4PVPlacement(rotBin,position,"DeltaE_BINA",deltaE_B_log[nrde], Mars_phs,false,NumBina[gdE]);  gdE++;
	       deltaE_B_phs[gdE]=new G4PVPlacement(rotBin, G4ThreeVector(xb, yb, zb),"DeltaE_BINA", deltaE_B_log[nrde], Mars_vac_phs,false,NumBina[gdE]);  gdE++; 
	     }
	   rotBin++;
	 }
     }
   if(ballvisible)
     {
       G4VisAttributes* kolorek= new G4VisAttributes(G4Colour(0.5,0.,0.5));
       Edet_B_log[0]  ->SetVisAttributes(kolorek);
       G4VisAttributes* kolorek1= new G4VisAttributes(G4Colour(0.5,1.0,0.5));
       Edet_B_log[1] ->SetVisAttributes(kolorek1);
       G4VisAttributes* kolorek2= new G4VisAttributes(G4Colour(0.,1.0,1.));
       Edet_B_log[2] ->SetVisAttributes(kolorek2);
       G4VisAttributes* kolorek3= new G4VisAttributes(G4Colour(1.,1.0,0.));
       Edet_B_log[3] ->SetVisAttributes(kolorek3);
       G4VisAttributes* kolorek4= new G4VisAttributes(G4Colour(0.5,1.0,1.));
       Edet_B_log[4] ->SetVisAttributes(kolorek4);
       G4VisAttributes* kolorek5= new G4VisAttributes(G4Colour(1.,0.5,0.75));
       Edet_B_log[5] ->SetVisAttributes(kolorek5);
       G4VisAttributes* kolorek7= new G4VisAttributes(G4Colour(1.,0.75,0.5));
       Edet_B_log[6] ->SetVisAttributes(kolorek7);
       G4VisAttributes* kolorek6= new G4VisAttributes(G4Colour(0.5,0.5,1.));
       G4VisAttributes* bilej= new G4VisAttributes(G4Colour(1.,1.,1.));
       pipe_log->SetVisAttributes(bilej);
       for (i=0;i<7;i++)
	 {
	   deltaE_B_log[i] ->SetVisAttributes(kolorek6);;
	 }
     }
   else
     {
       for (i=0;i<7;i++)
	 {
	   Edet_B_log[i] ->SetVisAttributes (G4VisAttributes::Invisible);
	   deltaE_B_log[i]->SetVisAttributes (G4VisAttributes::Invisible);
	 }
     }
 }  
 
  //--------------------------------
  //         Delta E NOWE
  //--------------------------------
if(deIs>0)
{
 int ide=0,nde ;
 //
  walec = new G4Tubs("walec",0.,deHoleOut,20.*cm,0.,2*M_PI);
  deltae_foil_sol=new G4Box("deltae_foil_sol",
		 dedimX,dedimY,dedimf);
  deltae_W_sol=new G4Box("deltae_W_sol",dedimX,dedimY,dedimZ);
 
  ///////liczenie ilosci pocietych scyntylatorow //////
  int cut;
  cut =int(deHoleOut/(2*(dedimX+deSepar))+1);

  for(i=nde_min;i<cut;i++)
  {
    // position liczone w celu odpowiedniego odciecia walca od detektora
   position = G4ThreeVector(-(2*i+1)*(dedimX+deSepar),0.,0.); 
   deltae_foilc_sol[i] = new G4SubtractionSolid("deltae_foilc_sol",
       deltae_foil_sol,walec,0,position);
   deltae_foil_log[i] = new G4LogicalVolume(deltae_foilc_sol[i], Al, 
      "deltae_foilc_sol");  
   deltae_Wc_sol[i] = new G4SubtractionSolid("deltae_Wc_sol",
       deltae_W_sol,walec,0,position);
   deltae_W_log[i] = new G4LogicalVolume(deltae_Wc_sol[i], BC400, 
       "deltae_wc_sol");  
  }
  deltae_foil_log[i] =new G4LogicalVolume(deltae_foil_sol,Al,"deltae_foil_log");
  deltae_W_log[i] =new G4LogicalVolume(deltae_W_sol,BC400,"deltae_W_log");

 
   if(deVisible)
   {
     G4VisAttributes* DeltaE= new G4VisAttributes(G4Colour(1.0,0.,0.));
     for(i=0;i<=cut;i++){  deltae_W_log[i]->SetVisAttributes(DeltaE);  }   }
   else
   { for(i=0;i<=cut;i++){ deltae_W_log[i]->SetVisAttributes(G4VisAttributes::Invisible); } }

   G4RotationMatrix* rot0;
   rot0 = new G4RotationMatrix(0,0,M_PI);
   ide =0;
   for (nde=0; nde < nde_max/2 ; nde++)
   {
     // right side
     // front foil
     position = G4ThreeVector(deXplace+((2*nde+1)*(dedimX+deSepar)), 
         deYplace, deZplace-dedimZ-dedimf-eps_separ);
     deltae_foil_phs[ide] = new G4PVPlacement(0,position,"DE_foil_phs",
           deltae_foil_log[(nde<cut)?nde:cut], Mars_phs,false,ide);
     // back foil
     position = G4ThreeVector(deXplace+((2*nde+1)*(dedimX+deSepar)), 
         deYplace, deZplace+dedimZ+dedimf+eps_separ);
     deltae_foil_phs[ide] = new G4PVPlacement(0,position,"DE_foil_phs",
           deltae_foil_log[(nde<cut)?nde:cut], Mars_phs,false,ide+nde_max/2);
     // scintillator
     position = G4ThreeVector(deXplace+((2*nde+1)*(dedimX+deSepar)), 
         deYplace, deZplace);
     deltae_W_phs[ide] = new G4PVPlacement(0,position,"DeltaE_W_phs",
          deltae_W_log[(nde<cut)?nde:cut], Mars_phs ,false, ide); ++ide;

     // Left side
     // front foil
     position = G4ThreeVector(-deXplace-((2*nde+1)*(dedimX+deSepar)),
            deYplace, deZplace-dedimZ-dedimf-eps_separ);
     deltae_foil_phs[ide] = new G4PVPlacement(rot0,position,"DE_foil_phs",
        deltae_foil_log[(nde<cut)?nde:cut], Mars_phs,false,ide);
     // back foil
     position = G4ThreeVector(-deXplace-((2*nde+1)*(dedimX+deSepar)),
            deYplace, deZplace+dedimZ+dedimf+eps_separ);
     deltae_foil_phs[ide] = new G4PVPlacement(rot0,position,"DE_foil_phs",
        deltae_foil_log[(nde<cut)?nde:cut], Mars_phs,false,ide+nde_max/2);
     // scintillator
     position = G4ThreeVector(-deXplace-((2*nde+1)*(dedimX+deSepar)),
            deYplace, deZplace);
     deltae_W_phs[ide] = new G4PVPlacement(rot0,position,"DeltaE_W_phs",
          deltae_W_log[(nde<cut)?nde:cut], Mars_phs,false, ide); ++ide;
   }

} 
 
  ///////////////////////
  ///  Degrader Foil
  ///////////////////////
 if(foilIs>0)
 {
  walec = new G4Tubs("walec",0.,foilholeOut,20.*cm,0.,2*M_PI);
  foil_sol=new G4Box("folia_sol",foildimX,foildimY,foildimZ);
  position = G4ThreeVector(0.,0.,0.); 
  foil_min_sol = new G4SubtractionSolid("folia_min_sol",foil_sol,walec,0,position);
  foil_log = new G4LogicalVolume(foil_min_sol, SSteel, "dfolia_log"); 
  position = G4ThreeVector(foilXplace,foilYplace,foilZplace);   
  foil_phs = new G4PVPlacement(0,position,"Foil_phs",foil_log,Mars_phs,false,0);
  if(foilVisible)
  {
    G4VisAttributes* kolorf= new G4VisAttributes(G4Colour(0.,1.0,1.));
    foil_log ->SetVisAttributes(kolorf);
  }
  else
  { foil_log->SetVisAttributes (G4VisAttributes::Invisible);  }
 }

  ////////////////////////////////////
  ///        Wall (E)
  ////////////////////////////////////

  int WallExist=1;
  int WnSide = 1;			//number of segments
  int WnZplanes = 2;			//number  of planes along z

  walec = new G4Tubs("walec",wallholeIn,wallholeOut,20.*cm,0.,2*M_PI);

  const G4double WZplanes_s[2] = {-wallDimcX , wallDimcX};	//z coordinate of  corners
  const G4double Wr_in_s[2] = {wallRIn , wallRIn};		//r_in coordinate of  corners
  const G4double Wr_out_s[2] = {wallROut , wallROut};	//r_out coordinate of  corners
  const G4double Wr_in_s1[2] = {wallRIn-dedimf , wallRIn-dedimf};
  const G4double Wr_out_s1[2] = {wallROut+dedimf , wallROut+dedimf};
  
  wallf_sol = new G4Polyhedra("Wallf_solid",wallAng0,wallAng,WnSide,WnZplanes,WZplanes_s,Wr_in_s1,Wr_out_s1);
 for (i=0;i<10;i++)  wallf_log[i] = new G4LogicalVolume(wallf_sol, Al, "Wallf_logical");

  wall_sol = new G4Polyhedra("Wall_solid",wallAng0,wallAng,WnSide,WnZplanes,WZplanes_s,Wr_in_s,Wr_out_s);
 for (i=0;i<10;i++) wall_log[i] = new G4LogicalVolume(wall_sol, BC408, "Wall_logical");

  G4RotationMatrix* rotFin = new G4RotationMatrix();
  rotFin->rotateY(M_PI/2.);
  position = G4ThreeVector((wallROut+wallRIn)/2.,0.,0.);
  //position = G4ThreeVector(wallROut,0.,0.);  

  //foil
  wallf_C_sol=new G4SubtractionSolid("Obc_specjalne3",wallf_sol,walec,rotFin,position);
  for (i=0;i<2;i++)  wallf_C_log[i] = new G4LogicalVolume(wallf_C_sol, Al, "Wallf_Center_logical");

  wallf_L_sol = new G4Box("Wallf_L_solid",wallLdimX+dedimf,wallLdimY+dedimf,wallLdimZ+dedimf);
  for(i=0;i<10;i++)  wallf_L_log[i] = new G4LogicalVolume(wallf_L_sol , Al,"Wallf_L_log");
  //scintillator 
  wall_C_sol=new G4SubtractionSolid("Obciecie_specjalne3",wall_sol,walec,rotFin,position);
  for (i=0;i<2;i++)  wall_C_log[i] = new G4LogicalVolume(wall_C_sol, BC408, "Wall_Center_logical");

  wall_L_sol = new G4Box("Wall_L_solid",wallLdimX,wallLdimY,wallLdimZ);
  for(i=0;i<10;i++) wall_L_log[i] = new G4LogicalVolume(wall_L_sol , BC408,"Wall_L_log");

 G4double ly, lY, lz, lZ, specdimy, specdimY;
 G4double kAng, kSepar, fi, minfi;
 G4double placeSpecZ, placeSpecY;
 
 specdimY =14.*cm;
 
 kAng = wallAng - wallAng0;
 kSepar = atan2(2.*wallSepar,wallROut);
 ang = kAng + kSepar;
 fi = wallNCenter * ang ;
 minfi = (90.*deg)-fi;
 ly = wallRIn * cos(minfi) + wallSepar;
 lY = wallROut * cos(minfi) + wallSepar;
 lz = wallRIn * sin(minfi);
 lZ = wallROut * sin(minfi);
 placeSpecZ = (lZ - lz)/2. + lz;
 specdimy = walldimS - (lY - ly) ;
 placeSpecY = specdimY/2.+(specdimY-specdimy)/4. +ly + wallYplace;

  wallf_Spec_sol = new G4Trap("Wallf_Spec_solid",2.*(wallLdimX+dedimf),2.*(wallLdimZ+dedimf),(walldimS+dedimf),specdimy+dedimf);
for (i=0;i<2;i++)  wallf_Spec_log[i] = new G4LogicalVolume(wallf_Spec_sol, Al, "Wall_Spec_logical"); 
  wall_Spec_sol = new G4Trap("Wall_Spec_solid",2.*wallLdimX,2.*wallLdimZ,walldimS,specdimy);
for (i=0;i<2;i++)  wall_Spec_log[i] = new G4LogicalVolume(wall_Spec_sol, BC408, "Wall_Spec_logical"); 

  G4double zs;
  zs = placeSpecY +specdimY/2.-(specdimY-specdimy)/4. + wallLdimY + 2.*wallSepar;
  //------------------------------
  //		Wall
  //------------------------------
  G4RotationMatrix* rotDo0= new G4RotationMatrix[11];
  G4RotationMatrix* rotUp0 = new G4RotationMatrix[11];
  G4RotationMatrix* rotDo5 = new G4RotationMatrix[11];
  G4RotationMatrix* rotUp5 = new G4RotationMatrix[11];
  G4RotationMatrix* rot5p = new G4RotationMatrix[11];
  rot5p = new G4RotationMatrix(0,0,0);
  ang_dif=ang;
  ang = 0.;
double r_ang;
int w=200;

 if (WallExist)
 {
 G4ThreeVector positions=G4ThreeVector(0,0,0);
// G4ThreeVector position;
  for (i=0;i<10;i++)	
  {
    r_ang=ang;// *M_PI/180.;

    if (i == 0)
      { //upper part
	rotDo0=new G4RotationMatrix(-M_PI/2.,-M_PI/2.,-M_PI/2.+r_ang);
	position = G4ThreeVector(wallXplace, -wallSepar,wallZplace); //foil
	wallf_phs[2*i] = new G4PVPlacement(rotDo0,position,"EFoil_Wall_phs0" , wallf_C_log[0] , Mars_phs,false,w);
	wall_phs[2*i] = new G4PVPlacement(rot5p,positions,"EDet_Wall_phs0" , wall_C_log[0] , wallf_phs[2*i] ,false,w++);
	//gorne
	rotUp0=new G4RotationMatrix(-M_PI/2.,M_PI/2.,M_PI/2.+r_ang);
	position = G4ThreeVector(wallXplace, wallSepar,wallZplace); //foil
       	wallf_phs[2*i+1] = new G4PVPlacement(rotUp0,position,"EFoil_Wall_phs1" , wallf_C_log[1] , Mars_phs,false,w);
	wall_phs[2*i+1] = new G4PVPlacement(rot5p,positions,"EDet_Wall_phs1" , wall_C_log[1] , wallf_phs[2*i+1],false,w++);
      }
    
    
    if (i>0&&i<wallNCenter)//(i >0&&i<wallNCenter)
      {
	rot1[i-1] = new G4RotationMatrix(-M_PI/2.,-M_PI/2.,-M_PI/2.+r_ang);
	position = G4ThreeVector(wallXplace,wallYplace-wallSepar,wallZplace);
	wallf_phs[2*i] = new G4PVPlacement(rot1[i-1],position,"EFoil_Wall_phs0-5",wallf_log[2*i-2],Mars_phs,false,w);
	wall_phs[2*i] = new G4PVPlacement(rot5p,positions,"EDet_Wall_phs0-5",wall_log[2*i-2],wallf_phs[2*i],false,w++);
	
	rot2[i-1] = new G4RotationMatrix(-M_PI/2.,M_PI/2.,M_PI/2.+r_ang);
	position = G4ThreeVector(wallXplace,wallYplace+wallSepar,wallZplace);
	wallf_phs[2*i+1] = new G4PVPlacement(rot2[i-1],position, "EFoil_Wall_phs0-5",wallf_log[2*i-1],Mars_phs,false,w);
	wall_phs[2*i+1] = new G4PVPlacement(rot5p,positions, "EDet_Wall_phs0-5",wall_log[2*i-1],wallf_phs[2*i+1],false,w++);
	}
	

    if (i==wallNCenter)
      {  //trapezoid
	rotDo5 = new G4RotationMatrix(-M_PI/2.,-M_PI/2.,M_PI);
	position = G4ThreeVector(wallXplace,-placeSpecY, placeSpecZ);//47.38
	wallf_phs[2*i] = new G4PVPlacement(rotDo5,position,"EFoil_Wall_phs5",wallf_Spec_log[0] , Mars_phs,false,w);
	wall_phs[2*i] = new G4PVPlacement(rot5p,positions,"EDet_Wall_phs5",wall_Spec_log[0] ,wallf_phs[2*i] ,false,w++);
	rotUp5= new G4RotationMatrix(-M_PI/2.,M_PI/2.,0.);
	position = G4ThreeVector(wallXplace,placeSpecY, placeSpecZ);
	wallf_phs[2*i+1] = new G4PVPlacement(rotUp5,position,"EFoil_Wall_phs5",wallf_Spec_log[1] , Mars_phs,false,w);
	wall_phs[2*i+1] = new G4PVPlacement(rot5p,positions,"EDet_Wall_phs5",wall_Spec_log[1] ,wallf_phs[2*i+1] ,false,w++);
      }
    if (i>wallNCenter)
      {
	if (i>6) {zs +=2.*(wallLdimY + wallSepar);}
	position = G4ThreeVector(wallXplace, (-zs), placeSpecZ);
	wallf_phs[2*i] = new G4PVPlacement(rot5p,position,"EFoil_Wall_phs6",wallf_L_log[2*i-12] , Mars_phs,false,w);
	wall_phs[2*i] = new G4PVPlacement(rot5p,positions,"EDet_Wall_phs6",wall_L_log[2*i-12] ,wallf_phs[2*i],false,w++);
	position = G4ThreeVector(wallXplace, (zs), placeSpecZ);//72.317
	wallf_phs[2*i+1] = new G4PVPlacement(rot5p,position,"EFoil_Wall_phs6",wallf_L_log[2*i-11] , Mars_phs,false,w);
	wall_phs[2*i+1] = new G4PVPlacement(rot5p,positions,"EDet_Wall_phs6",wall_L_log[2*i-11] , wallf_phs[2*i+1],false,w++);
      }
    ang += ang_dif;  	//rotate angle
  }
 }
 Seen();
 return Mars_phs;
}
  /*
void Bina_DetectorConstruction::SetParamUpdate() 
  {G4RunManager::GetRunManager()->SetUserAction(new Bina_PrimaryGeneratorAction(this)); }
*/