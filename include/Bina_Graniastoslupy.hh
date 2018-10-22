#ifndef G4GRANIASTOSLUPY_HH
#define G4GRANIASTOSLUPY_HH

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"///
#include "G4Polyhedra.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"

class polozenia
{
	public:
	G4int typ;
	G4double x,y,z;
	G4double rx,ry,rz;
	G4double r;
	polozenia();

};

class wymiary
{
  public:

  private:
  G4double a1,a2,a1m,a2m;             //wymiary podstawy gornej , dolnej i dolnej cienkiego scyntylatora
  G4double h1,h2,h1m,h2m;             //wysokosci podstaw
  G4double h,hm,r,rm;                      // wysokosci graniastoslupow i orleglosci od srodka
  G4double ak,bk,akr,bkr,gk,gkr;            //katy i radiany
  G4double bkm,akrm,bkrm,gkm,gkrm; //katy i radiany dla malego graniastoslupa
  G4double d, d1, q, c,b,zt;               //zmienne pomoocnicze
  G4double dm,d1m,qm,cm,bm,ztm;
  G4double dv1,dv2;                        // ciecia w graniastoslupach specjalnych 1,2,3
  G4double dv1m,dv2m;


 public:
  wymiary();
  wymiary(G4double ,G4double ,G4double ,G4double ,G4double );
  void licz();
  G4double get_a1() {return a1;};
  G4double get_a2() {return a2;};
  G4double get_a1m() {return a1m;};
  G4double get_a2m() {return a2m;};
  G4double get_h1() {return h1;};
  G4double get_h2() {return h2;};
  G4double get_h1m() {return h1m;};
  G4double get_h2m() {return h2m;};
  G4double get_h() {return h;};
  G4double get_hm() {return hm;};
  G4double get_r() {return r;};
  G4double get_rm() {return rm;};
  G4double get_ak() {return ak;};
  G4double get_bk() {return bk;};
  G4double get_akr() {return akr;};
  G4double get_bkr() {return bkr;};
  G4double get_gk() {return gk;};
  G4double get_gkr() {return gkr;};
  G4double get_d() {return d;};
  G4double get_d1() {return d1;};
  G4double get_q() {return q;};
  G4double get_c() {return c;};
  G4double get_b() {return b;};
  G4double get_zt() {return zt;};
  G4double get_dm() {return dm;};
  G4double get_d1m() {return d1m;};
  G4double get_qm() {return qm;};
  G4double get_cm() {return cm;};
  G4double get_bm() {return bm;};
  G4double get_ztm() {return ztm;};
  G4double get_dv1() {return dv1;};
  G4double get_dv2() {return dv2;};
  G4double get_dv1m() {return dv1m;};
  G4double get_dv2m() {return dv2m;};
//  void set();
};

G4VSolid*  przecinka(G4VSolid *,wymiary &,int ,int typ=0);
void planowanie(polozenia*,G4double R[],int);

// wys --- detektor gruby E ->1 ,cienki de ->2
// typ okresla specjalne ciecia
#endif
