Include "Angle.PRO"

Group {
  /* Defines Physical regions: 
  Air - Air gap
  Core - stator/rotor iron
  Surface_Shaft - Shaft Boundry
  Surface_Air - Air boundry with stator
  Surface_PN - Surface on P side at position N
  Surf_EdgeP - Surface on the P side
  Surface_AGC - Surface over the airgap
  */
  AirCore     = Region[ 10000 ];
  Core	  = Region[ 10004 ];
  AirT    = Region[ 10006 ];
  AirB    = Region[ 10002 ];
  SlotLiner = Region[ 10007 ];
  A1 = Region[ 10008 ];
  A2 = Region[ 10009 ];
  A3 = Region[ 10010 ];
  A4 = Region[ 10011 ];
  A5 = Region[ 10012 ];
  A6 = Region[ 10013 ];
  A7 = Region[ 10014 ];
  A8 = Region[ 10015 ];
  B1 = Region[ 10016 ];
  B2 = Region[ 10017 ];
  B3 = Region[ 10018 ];
  B4 = Region[ 10019 ];
  B5 = Region[ 10020 ];
  B6 = Region[ 10021 ];
  B7 = Region[ 10022 ];
  B8 = Region[ 10023 ];
  C1 = Region[ 10024 ];
  C2 = Region[ 10025 ];
  C3 = Region[ 10026 ];
  C4 = Region[ 10027 ];
  C5 = Region[ 10028 ];
  C6 = Region[ 10029 ];
  C7 = Region[ 10030 ];
  C8 = Region[ 10031 ];
  A_area = Region[{A1,A2,A3,A4,A5,A6,A7,A8}];
  B_area = Region[{B1,B2,B3,B4,B5,B6,B7,B8}];
  C_area = Region[{C1,C2,C3,C4,C5,C6,C7,C8}];
  Conductors = Region[{A_area,B_area,C_area}];
  SatorCore    = Region[ 10005 ];

  For i In {0:3}
    Surface_AGC~{i+1} = Region[ {(20000+(i))} ];
  	RotorBoundry += Region[ Surface_AGC~{i+1} ];
  EndFor

  For i In {0:1}
    Surface_R~{i} = Region[ {(20008+(i))} ];
    Surface_L~{i} = Region[ {(20012+(i))} ];
  	Sur_Periodic += Region[ {Surface_L~{i},Surface_R~{i}} ];
  EndFor
  Surface_R_2 = Region[ {(20027)} ];
  Surface_R_3 = Region[ {(20028)} ];
  Surface_L_2 = Region[ {(20029)} ];
  Surface_L_3 = Region[ {(20030)} ];
  Sur_Periodic += Region[ {Surface_L_2,Surface_R_2,Surface_L_3,Surface_R_3} ];
  
  Surface_Shaft = Region[ 20016 ]; //rotor shaft
  Dirchlet_1 = Region[ 20026 ];
  Surface_AGCB = Region[ 20019 ]; //airgap

  Periodic_Boundry_Air_R = Region[ 20020 ];
  Periodic_Boundry_Stator_Lower_R = Region[ 20021 ];
  Periodic_Boundry_Stator_Upper_R = Region[ 20022 ];
  Periodic_Boundry_Air_L = Region[ 20023 ];
  Periodic_Boundry_Stator_Lower_L = Region[ 20024 ];
  Periodic_Boundry_Stator_Upper_L = Region[ 20025 ];
  Sur_Periodic += Region[ {20020,20021,20022,20023,20024,20025} ];

  
  RotationCenter = Region[ 30002 ];
  MovingBand_PhysicalNb = Region[0];

  /* Abstract regions :
  */
  Sur_Dir = Region[ { Surface_Shaft,Dirchlet_1 } ];
  Surf_EdgeL = Region[ {Surface_L_0} ];
  CoreAll = Region[ {Core,AirCore} ];
  Surf_EdgeR = Region[ {Surface_R_0} ];  
  Surf_EdgeSL = Region[ {Periodic_Boundry_Stator_Lower_L,Periodic_Boundry_Stator_Upper_L} ];
  Surf_EdgeSR = Region[ {Periodic_Boundry_Stator_Lower_R,Periodic_Boundry_Stator_Upper_R} ];
  MB = MovingBand2D[MovingBand_PhysicalNb, Surface_AGCB, RotorBoundry, 4] ;
  Air     = Region[ {AirCore, AirB,MB,AirT} ];
  Iron     = Region[ {Core, SatorCore} ];//
  AirSurf = Region[ {} ];
  Volume     = Region[ {Air, Iron,SlotLiner,Conductors} ];
  Volume_Other     = Region[ {Air, SlotLiner,Conductors} ];
  VolumeOther     = Region[ {Air,AirSurf,Sur_Dir,Sur_Periodic,Surface_AGC_1,Surface_AGCB} ];
  VolumeAll	 = Region[ {Volume,AirSurf,Sur_Dir,Sur_Periodic,MB,Surface_AGC_1} ];
  RotorBoundryAux = Region[ {RotorBoundry,-Surface_AGC~{1}} ];
  Rotor_Moving = Region[ {Core,AirCore,AirB,RotorBoundryAux} ];
  AirBoundry = Region[{AirT,AirB,MB}];
  Dom_H_u_Mec = Region[{Core,AirCore,Surface_Shaft,Surface_L_0,Surface_R_0}];
  
  Right = Region[{Surf_EdgeSR,Surf_EdgeR,Surface_R_2,Surface_R_3,Periodic_Boundry_Air_R}];//,Surface_R_2,Surface_R_3
  Left = Region[{Surf_EdgeSL,Surf_EdgeL,Surface_L_2,Surface_L_3,Periodic_Boundry_Air_L}];//,Surface_L_1,Surface_L_2,Surface_L_3
}

Function {
  mu0 = 4.e-7 * Pi;
  murCore = 1000;
  Young = 200e9;
  Poisson = 0.3;
  RPM = 3000;
  Density = 7850;
  FE_Order = 2;
  Angle_Shift[] = -30-($Time*180/36);

  /*aplying material properties*/
  /*CURRENT LOAD*/
  A[] = current_multiplier*50*Cos[(240+Angle_Shift[])*(Pi/180)]; //43
  B[] = -current_multiplier*50*Cos[(120+Angle_Shift[])*(Pi/180)];
  C[] = current_multiplier*50*Cos[(0+Angle_Shift[])*(Pi/180)];
  
  
  Current[Region[{A_area}]] = A[]/SurfaceArea[];
  Current[Region[{B_area}]] = B[]/SurfaceArea[];
  Current[Region[{C_area}]] = C[]/SurfaceArea[];

  /*MAGNETIC PERMIABILITY*/
  nu [ Region[{Air}] ]  			= 1. / mu0;
  nu [ Region[{Surface_Shaft}] ] 	= 1 / ( mu0);
  nu [ Region[{RotorBoundry}] ] 	= 1. / (mu0*1000);
  nu [ Region[{Surface_AGCB}] ] 	= 1. / mu0;
  nu [ Region[{Periodic_Boundry_Air_L}] ]  	= 1. / (mu0);
  nu [ Region[{Periodic_Boundry_Air_R}] ]  	= 1. / (mu0);
  nu [ Region[{SlotLiner}] ]  	= 1. / (mu0);
  nu [ Region[{Conductors}] ]  	= 1. / (mu0);
  nu [ Region[{Surface_L_0}] ]  	= 1. / (mu0*1000);
  nu [ Region[{Surface_R_0}] ]  	= 1. / (mu0*1000);

  /*YOUNGS MODULOUS*/
  E [ Region[{Air}] ]  				= 1;
  E [ Region[{RotorBoundry}] ]  	= 1;
  E [ Region[{Surface_AGCB}] ]  	= 1;
  E [ Region[{Surf_EdgeL}] ]  		= 1;
  E [ Region[{Surf_EdgeR}] ]  		= 1;
  E [ Region[{Surface_Shaft}] ]     = 1;
  E [ Region[{Core}] ]     			= Young;

  /*POISSION NUMBER*/
  mu [ Region[{Air}] ]  			= 0.003;
  mu [ Region[{RotorBoundry}] ]  	= 0.003;
  mu [ Region[{Surface_AGCB}] ]  	= 0.003;
  mu [ Region[{Surf_EdgeL}] ]  		= 0.003;
  mu [ Region[{Surf_EdgeR}] ]  		= 0.003;
  mu [ Region[{Surface_Shaft}] ]    = 0.003;
  mu [ Region[{Core}] ]     		= Poisson; 

  /*Density*/
  DensityEval[Region[{Core}]] = 1;
  DensityEval[Region[{AirCore}]] = 0;
  Density[Region[{Core}]] = 7850;
  Density[Region[{AirCore}]] = 1.293;

    /*MAGNETIC PERMIABILITY*/
  Dam = 0.95;
  a1 = 1.35E-03; a2 = 1.378426182; a3 = 6.87E-06 ; a4 = 0.195520764;
  Bl[]     = (a1*($1^(a2)))/(1+a3*($1^a4));
  damage[] = (1-Dam*Exp[-Max[Bl[Min[$1,500]],0.01]*$2*1e-0]);

  Mat_core_b = {
      0, 0.030000,0.060000,0.120000,0.190000,0.230000,0.580000,
      0.770000,0.950000,1.130000,1.310000,1.490000,1.690000,
	  1.830000,1.920000,1.980000,2.030000,2.080000,2.150000,
	  2.210000,2.290000,2.330000};
  Mat_core_h = {
      0.000000,93.000000,165.000000,284.000000,399.000000,
	  457.000000,1104.000000,1594.000000,2306.000000,3606.000000,
	  6468.000000,12904.000000,26799.000000,49770.000000,
	  74770.000000,99770.000000,124770.000000,149770.000000,
	  189770.000000,229770.000000,279770.000000,304770.000000};
  Mat_core_b2 = Mat_core_b()^2;
  Mat_core_h2 = Mat_core_h()^2;
  Mat_core_nu = Mat_core_h() / Mat_core_b();
  Mat_core_nu(0) = Mat_core_nu(1);
  Mat_core_nu_b2  = ListAlt[Mat_core_b2(), Mat_core_nu()] ;
  nu_core[] = InterpolationAkima[ SquNorm[$1] ]{ Mat_core_nu_b2() }; // *(1/InterpolationAkima[Norm[$1]+1]{Dens()}) ;
  dnudb2_core[] = dInterpolationAkima[SquNorm[$1]]{ Mat_core_nu_b2() } ;
  h_core[] = nu_core[$1] * $1 ;
  h2_core[] = dnudb2_core[$1] * $1 ;
  dhdb_core[] = TensorDiag[1,1,1] * nu_core[$1#1] +
    2*dnudb2_core[#1] * SquDyadicProduct[#1]  ;
  /* Stator iron assumed nonlinear no damage */
  nu[ Region[{SatorCore}] ]  	= nu_core[$1];
  dhdb[ Region[{SatorCore}] ] 	= dhdb_core[$1];
  /* rotor iron assumed non linerar + damage */
  nu[ Region[{Core}] ]  		= nu_core[$1]*1/damage[Norm[h_core[$1] ],$2];
  dhdb[ Region[{Core}] ] 		= dhdb_core[$1]*1/damage[Norm[h2_core[$1] ],$2];
  CS [ Region[{Iron}]] = 0;
  
  RotatePZ[] = Rotate[ Vector[$X,$Y,$Z], 0, 0, $1 ] ;
  TM[] = ( SquDyadicProduct[$1] - SquNorm[$1] * TensorDiag[0.5, 0.5, 0.5] ) / (1/nu[]) ;
  
  /*Angle calculations*/
  Angle[] = Atan2[Y[],X[]];
  
    /*Stress calculations*/ 
    /*Stress Constants*/ 
    a[] = E[]/(1.-mu[]^2);
    c[] = E[]*mu[]/(1.-mu[]^2);
    b[] = E[]/2./(1.+mu[]);
    /*2D Stress Tensor*/ 
    C_xx[] = Tensor[ a[],0  ,0  ,    0  ,b[],0  ,    0  ,0  ,b[] ];
    C_xy[] = Tensor[ 0  ,c[],0  ,    b[],0  ,0  ,    0  ,0  ,0   ];
    C_yx[] = Tensor[ 0  ,b[],0  ,    c[],0  ,0  ,    0  ,0  ,0   ];
    C_yy[] = Tensor[ b[],0  ,0  ,    0  ,a[],0  ,    0  ,0  ,b[] ];
	
    Mesies[] = Sqrt[((CompX[  C_xx[]*$1 + C_xy[]*$2 ])^2+(CompY [ C_yx[]*$1 + C_yy[]*$2 ])^2+(CompX [ C_yx[]*$1 + C_yy[]*$2 ])^2+(CompY [ C_yx[]*$2 + C_xx[]*$1])^2)];

    /*Force Calculation and conversion to pressure (done seperate for unkown reason)*/ 
    omega[] = RPM / 60;
    rad[] = Hypot[X[], Y[]];
    force_x[] = Density[] * omega[]^2 * rad[] * Cos[Angle[]];
    force_y[] = Density[] * omega[]^2 * rad[] * Sin[Angle[]];
}

/*CONSTRAINTS*/
Constraint {
  { Name Shaft ; Type Assign;
    Case {
      /*SHAFT*/
      { Region Sur_Dir ; Value 0.; }

      /*PERIODICITY*/
      { Region Right; SubRegion Region[{Sur_Dir, RotationCenter}]; Type Link;
          RegionRef Left; SubRegionRef Region[{Sur_Dir, RotationCenter}];
          Coefficient ((1%2)?-1:1) ;
	  Function RotatePZ[2*Pi/4]; }
	  
        { Region Right; Type Link; RegionRef Left;
          Coefficient (1%2)?-1:1 ;
	  Function RotatePZ[2*Pi/4]; }

    /*STRESS TENSOR*/
      For k In {1:4-1}
        { Region Surface_AGC~{k+1} ;
	    SubRegion Surface_AGC~{(k!=4-1)?k+2:1}; Type Link;
        RegionRef Surface_AGC_1; SubRegionRef Surface_AGC_2;
        Coefficient ((1%2)?-1:1)^(k);
	    Function RotatePZ[-k*1*2*Pi/4]; }
      EndFor
         }
  }
  { Name ForceConstraint;
    Case{
      { Region Surface_AGC_1 ; Value 1. ; }
    }
  }
  /*
  deflection in x BCs
  sets shft deflection in x to 0
  */ 
  { Name Displacement_x;
  Case {
        {Region Sur_Dir ; Type Assign ; Value 0; }
    }
  }
  
  /*
  deflection in y BCs
  sets shft deflection in y to 0
  */ 
  { Name Displacement_y;
    Case {
        {Region Sur_Dir ; Type Assign ; Value 0; }
    }
  } 
  /*no distance constraints*/
  { Name PassOver ; Type Assign;
    Case {
         }
  }
  }

FunctionSpace {
  { Name Hgrad_v; Type Form0;
    BasisFunction {
      { Name se; NameOfCoef ae; Function BF_Node;
        Support Region[{ VolumeAll, RotorBoundryAux }]; Entity NodesOf[All]; }
    }
    Constraint {
      { NameOfCoef ae; EntityType NodesOf;
        NameOfConstraint Shaft; }
    }	
  }
  { Name StressTensor ; Type Form0 ;
    BasisFunction {
      { Name sn ; NameOfCoef un ; Function BF_GroupOfNodes ;
        Support MB ; Entity GroupsOfNodesOf[ Surface_AGC_1 ] ; }
    }
    Constraint {
      { NameOfCoef un ; EntityType GroupsOfNodesOf ; NameOfConstraint ForceConstraint; }
    }
  }
  /*
  Mechnaical function space in x with Displacement_x BCs
  Form 2E - divergace field
  Support - mechanical domain
  constraint - displacement in x
  */ 
{ Name H_ux_Mec ; Type Form0 ;
    BasisFunction {
      { Name sxn ; NameOfCoef uxn ; Function BF_Node ;
        Support Dom_H_u_Mec ; Entity NodesOf[ All ] ; }
     If ( FE_Order == 2 )
        { Name sxn2 ; NameOfCoef uxn2 ; Function BF_Node_2E ;
          Support Dom_H_u_Mec; Entity EdgesOf[ All ] ; }
     EndIf
    }
    Constraint {
      { NameOfCoef uxn ;
        EntityType NodesOf ; NameOfConstraint Displacement_x ; }
      If ( FE_Order == 2 )
      { NameOfCoef uxn2 ;
        EntityType EdgesOf ; NameOfConstraint Displacement_x ; }
      EndIf
    }
  }
  /*
  Mechnaical function space in y with Displacement_y BCs
  Form 2E - divergace field
  Support - mechanical domain
  constraint - displacement in y
  */ 
  { Name H_uy_Mec ; Type Form0 ;
    BasisFunction {
      { Name syn ; NameOfCoef uyn ; Function BF_Node ;
        Support Dom_H_u_Mec ; Entity NodesOf[ All ] ; }
     If ( FE_Order == 2 )
        { Name syn2 ; NameOfCoef uyn2 ; Function BF_Node_2E ;
          Support Dom_H_u_Mec; Entity EdgesOf[ All ] ; }
     EndIf
    }
    Constraint {
      { NameOfCoef uyn ;
        EntityType NodesOf ; NameOfConstraint Displacement_y ; }
      If ( FE_Order == 2 )
      { NameOfCoef uyn2 ;
        EntityType EdgesOf ; NameOfConstraint Displacement_y ; }
      EndIf
    }
  }
  /*
  Mechnaical lagrange function implementing periodic boundry conitions 
  assuming field normal for boundry
  Form 0 - scaler field
  Support - periodic_R1/periodic_L1
  constraint - 
  */ 
    { Name Lagrange ; Type Form0 ;
    BasisFunction {
      { Name lag_sn ; NameOfCoef lag_un ; Function BF_Node ;
        Support Surface_R_0 ; Entity NodesOf[ All ] ; }
    }
    Constraint {
    }
  }
    { Name LagrangeM ; Type Form0 ;
    BasisFunction {
      { Name lagm_sn ; NameOfCoef lagm_un ; Function BF_Node ;
        Support Surface_L_0 ; Entity NodesOf[ All ] ; }
    }
    Constraint {
    }
  }
  /*
  Mechnaical function space for VonMisies stress with Sig BCs
  Form 0 - scaler field
  Support - VolumeMech
  constraint - 
  */ 
  { Name H_Sig_Mec; Type Form0;
    BasisFunction {
      { Name sSign; NameOfCoef uSign; Function BF_Node;
      Support Core; Entity NodesOf[ All ]; }
    }
    Constraint {
    }
  }
  /*Distance domain*/
  { Name Dist_v; Type Form0;
    BasisFunction {
      { Name sise; NameOfCoef dise; Function BF_Node;
        Support Region[{ VolumeAll }]; Entity NodesOf[All]; }
    }
    Constraint {
      { NameOfCoef dise; EntityType NodesOf;
        NameOfConstraint PassOver; }
    }	
  }
}


Jacobian {
  { Name JacVol ;
    Case { 
           { Region All ; Jacobian Vol ; }
    }
  }
  { Name Sur ;
    Case { 
           { Region All ; Jacobian Sur ; }
    }
  }
    { Name JLin ;
    Case {
      { Region All ; Jacobian Lin ; }
    }
  }
}

Integration {
  { Name Int ;
    Case { { Type Gauss ;
             Case { 
				   { GeoElement Point   ; NumberOfPoints 1; }
				   { GeoElement Line   ; NumberOfPoints 4; }
				   { GeoElement Triangle   ; NumberOfPoints 4; }
                   { GeoElement Quadrangle ; NumberOfPoints 4; }
                   { GeoElement Tetrahedron; NumberOfPoints 4; }
                   { GeoElement Hexahedron ; NumberOfPoints 6; }
                   { GeoElement Prism      ; NumberOfPoints 9; }
	}
      }
    }
  }
}
//needds work
Formulation{
	{ Name Magnetics; Type FemEquation;
		Quantity{
		{Name a; Type Local; NameOfSpace Hgrad_v;}
		{Name un ; Type Local ; NameOfSpace StressTensor ; }
	    {Name dis; Type Local; NameOfSpace Dist_v;}
		}
	
	Equation{
	    Integral{[  0*Dof{d a} , {d a} ];
				In RotorBoundryAux; Jacobian Sur; Integration Int;}
	    Integral{[CS[],{a}];
				In AirSurf; Jacobian JacVol; Integration Int;}
		Integral{[ nu[]*Dof{d a},{d a}];
				In Air; Jacobian JacVol; Integration Int;}
		Integral{[ nu[]*Dof{d a},{d a}];
				In SlotLiner; Jacobian JacVol; Integration Int;}
		Integral{ [ nu[]*Dof{d a} ,{d a}];
				In Conductors; Jacobian JacVol; Integration Int;}       // flux potential from point sorce
	    Integral{ [ Current[] , {a} ];
				In Conductors; Jacobian JacVol; Integration Int;}       // flux potential from point sorce
	   
		Integral{[ nu[{d a},{dis}]*Dof{d a} ,{d a}];
				In Core; Jacobian JacVol; Integration Int;}
		Integral{[ dhdb[{d a},{dis}]*Dof{d a} ,{d a}];
				In Core; Jacobian JacVol; Integration Int;}
		Integral{[ - dhdb[{d a},{dis}]*{d a} ,{d a}];
				In Core; Jacobian JacVol; Integration Int;}
		Integral{[ nu[{d a}]*Dof{d a} ,{d a}];
				In SatorCore; Jacobian JacVol; Integration Int;}
		Integral{[ dhdb[{d a}]*Dof{d a} ,{d a}];
				In SatorCore; Jacobian JacVol; Integration Int;}
		Integral{[ - dhdb[{d a}]*{d a} ,{d a}];
				In SatorCore; Jacobian JacVol; Integration Int;}
		Integral{[ 0 * Dof{un} , {un} ];
				In VolumeAll ; Jacobian JacVol; Integration Int;}
	}
}
  /*Distance weak formulation*/
  { Name DistanceCalcuations; Type FemEquation;
    Quantity{
	   {Name dis; Type Local; NameOfSpace Dist_v;}
	}
	Equation{
		Integral{[-Field[XYZ[]],{dis}];
			In VolumeAll; Jacobian JacVol; Integration Int;}
		Integral{[Dof{dis},{dis}];
			In VolumeAll; Jacobian JacVol; Integration Int;}
	}
  }
/*Defenition of Deflection Problem*/ 
  { Name Elast_u; Type FemEquation;
    Quantity{
      { Name ux  ; Type Local ; NameOfSpace H_ux_Mec ; }   			// Deflection in x (radial)
      { Name uy  ; Type Local ; NameOfSpace H_uy_Mec ; }   			// Deflection in y (rotational)
	  { Name lag_lambda  ; Type Local ; NameOfSpace Lagrange ; }	// LAGRANGE MULTILIER LEFT SIDE
	  { Name lag_lambdaM  ; Type Local ; NameOfSpace LagrangeM ; }  // LAGRANGE MULTILIER RIGHT SIDE
	}
    Equation {
	  // Stress tensors ------------------------------------------------------------
      Integral { [ -C_xx[] * Dof{d ux}, {d ux} ] ;
        In CoreAll ; Jacobian JacVol ; Integration Int ; }
      Integral { [ -C_xy[] * Dof{d uy}, {d ux} ] ;
        In CoreAll ; Jacobian JacVol ; Integration Int ; }
      Integral { [ -C_yx[] * Dof{d ux}, {d uy} ] ;
        In CoreAll ; Jacobian JacVol ; Integration Int ; }
      Integral { [ -C_yy[] * Dof{d uy}, {d uy} ] ;
        In CoreAll ; Jacobian JacVol ; Integration Int ; }		

      //// Lagrange multiplier------------------------------------------------------
      // int (lag_lambda * (n * test(u'))) 
	  // where: lag_lambda is the Dof, n is the normal (Cos(angle normal to slave)) to slave and 
	  // u' the test function
      Integral { [ - Cos[Angle[]-Pi/2]* Dof{lag_lambda},  {ux} ];
        In Surface_R_0 ; Jacobian Sur ; Integration Int ; }
      Integral { [ - Sin[Angle[]-Pi/2]* Dof{lag_lambda},  {uy} ];
        In Surface_R_0 ; Jacobian Sur ; Integration Int ; }
      Integral { [ - Cos[Angle[]+Pi/2]* Dof{lag_lambdaM},  {ux} ];
        In Surface_L_0 ; Jacobian Sur ; Integration Int ; }
      Integral { [ - Sin[Angle[]+Pi/2]* Dof{lag_lambdaM},  {uy} ];
        In Surface_L_0 ; Jacobian Sur ; Integration Int ; }

	  // int ((n * u) * test(lag_lambda)) 
	  // where: , n is the normal to slave and u the Dof
	  // test(lag_lambda) is test function of lag_lambda
      Integral { [ -Dof{ux} * Cos[Angle[]-Pi/2]  , {lag_lambda} ];
        In Surface_R_0 ; Jacobian Sur ; Integration Int ; }
      Integral { [ -Dof{uy} * Sin[Angle[]-Pi/2]  , {lag_lambda} ];
        In Surface_R_0 ; Jacobian Sur ; Integration Int ; }
      Integral { [ -Dof{ux} * Cos[Angle[]+Pi/2]  , {lag_lambdaM} ];
        In Surface_L_0 ; Jacobian Sur ; Integration Int ; }
      Integral { [ -Dof{uy} * Sin[Angle[]+Pi/2]  , {lag_lambdaM} ];
        In Surface_L_0 ; Jacobian Sur ; Integration Int ; }
		
	  // Forces -------------------------------------------------------------------
      Integral { [ force_x[] , {ux} ];
        In CoreAll ; Jacobian JacVol ; Integration Int ; }
      Integral { [ force_y[] , {uy} ];
        In CoreAll ; Jacobian JacVol ; Integration Int ; }
	}
  }
}

Resolution {
  { Name Magnetics;
    System {
      { Name Sys_Mag; NameOfFormulation Magnetics; }
      { Name Sys_Mec; NameOfFormulation Elast_u; }               	//Links to Elastic problem
	  { Name distanceimport; NameOfFormulation DistanceCalcuations; }	//Links to Distacne problem
    }
    Operation {
	  /*Solves Distance problem*/ 
	  GmshRead["Distance_mini.pos",1];
	  Generate[distanceimport];
	  Solve[distanceimport];
	  SetTime[0];
	  /*Solves stress problem*/ 
	  //If (Angle < 5)
	  /*Solves mag problem*/ 
	  If (Angle == 0)
		  InitMovingBand2D[MB];
		  MeshMovingBand2D[MB];
		  ChangeOfCoordinates[ NodesOf[Rotor_Moving] , RotatePZ[Pi/8]];//Surf_EdgeAL
		  MeshMovingBand2D[MB] ;
		  Generate[Sys_Mag]; Solve[Sys_Mag];
		  Generate[Sys_Mag]; GetResidual[Sys_Mag,$res0];
		  Evaluate[$res = $res0,$iteration = 0];
		  Print[{$iteration, $res, $res / $res0},Format "Residual %03g: abs %14.12e rel %14.12e"];
		  While[$res > 1e-6 && $res / $res0 > 1e-6 && $iteration < 10]{
			Solve[Sys_Mag]; Generate[Sys_Mag]; GetResidual[Sys_Mag, $res];
			Evaluate[ $iteration = $iteration + 1 ];
			Print[{$iteration, $res, $res / $res0},Format "Residual %03g: abs %14.12e rel %14.12e"];
		  }
		  PostOperation[Torque];
		  SetTime[0];
		  PostOperation[Torque_ripple];
		  Evaluate[$iteration_2 = 1];
		  While[$iteration_2 < 4]{
			  ChangeOfCoordinates[ NodesOf[Rotor_Moving] , RotatePZ[Pi/64]];//Surf_EdgeAL
			  SetTime[$iteration_2];
			  MeshMovingBand2D[MB] ;
			  Generate[Sys_Mag]; Solve[Sys_Mag];
			  Generate[Sys_Mag]; GetResidual[Sys_Mag,$res0];
			  Evaluate[$res = $res0,$iteration = 0];
			  Print[{$iteration, $res, $res / $res0},Format "Residual %03g: abs %14.12e rel %14.12e"];
			  While[$res > 1e-6 && $res / $res0 > 1e-6 && $iteration < 10]{
				Solve[Sys_Mag]; Generate[Sys_Mag]; GetResidual[Sys_Mag, $res];
				Evaluate[ $iteration = $iteration + 1 ];
				Print[{$iteration, $res, $res / $res0},Format "Residual %03g: abs %14.12e rel %14.12e"];
			  }
			Evaluate[$iteration_2 = $iteration_2 + 1];
		    PostOperation[Torque_ripple];
		  }
	  EndIf
	  
	  SetTime[0];
	  ChangeOfCoordinates[ NodesOf[Rotor_Moving] , RotatePZ[-3*Pi/74]];//Surf_EdgeAL
	  ChangeOfCoordinates[ NodesOf[Rotor_Moving] , RotatePZ[-Pi/8]];//Surf_EdgeAL
	  MeshMovingBand2D[MB] ;
	  Generate[Sys_Mag]; Solve[Sys_Mag];
	  Generate[Sys_Mag]; GetResidual[Sys_Mag,$res0];
	  Evaluate[$res = $res0,$iteration_2 = 0];
	  Print[{$iteration_2, $res, $res / $res0},Format "Residual %03g: abs %14.12e rel %14.12e"];
	  While[$res > 1e-6 && $res / $res0 > 1e-6  && $res / $res0 <= 1 && $iteration_2 < 10]{
		Solve[Sys_Mag]; Generate[Sys_Mag]; GetResidual[Sys_Mag, $res];
		Evaluate[ $iteration_2 = $iteration_2 + 1 ];
		Print[{$iteration_2, $res, $res / $res0},Format "Residual %03g: abs %14.12e rel %14.12e"];
	  }
	  PostOperation[Map];
	  SaveSolution[Sys_Mag];
	  //Else
	  Generate[Sys_Mec];
	  Solve[Sys_Mec];
	  PostOperation[pos];
	  SaveSolution[Sys_Mec];
	  Exit;
	  //EndIf
    }
  }
}

PostProcessing {
  { Name Sys_Mag_2D; NameOfFormulation Magnetics;
    Quantity {
      { Name a; Value {Term {[ {a} ]; In VolumeAll; Jacobian JacVol; }} }
      { Name b; Value {Term {[Norm[ {d a} ]];In VolumeAll; Jacobian JacVol; }} } 
      { Name bLine; Value {Term {[ Norm[{d a}] ];In VolumeAll; Jacobian JLin; }} } 
      { Name Den; Value {Term {[DensityEval[ElementNum[]] ];In CoreAll; Jacobian JacVol; }} }
	  { Name AREA; Value {Term {[SurfaceArea[] ];In VolumeAll; Jacobian JacVol; }} }
	  { Name un ; Value { Local { [ {un} ] ; In VolumeAll; Jacobian JacVol ; } } }
	  { Name CurrentStream ; Value { Local { [ Current[]/10^6 ] ; In Conductors; Jacobian JacVol ; } } }
      { Name f ; Value { Integral { [ - TM[{d a}] * {d un} ] ;
              In MB ; Jacobian JacVol ; Integration Int ; } } }
	  { Name TorqueMaxwell ; Value { Integral {[CompZ[ -XYZ[] /\ (TM[{d a}] * XYZ[])* 2*Pi*0.05/SurfaceArea[]]];
              In Volume ; Jacobian JacVol ; Integration Int ; } } }    
	{ Name Torque_vw ; Value {Integral { Type Global ;
           [ CompZ[ 0.5 * nu[] * XYZ[] /\ VirtualWork[{d a}] ] * 0.1 ];
           In ElementsOf[Surface_AGC_1, OnOneSideOf MB];
	   Jacobian JacVol ; Integration Int ; }
       }}
    }
}
/*Defenition of Mechanical equations*/ 
  { Name Elast_u_2D; NameOfFormulation Elast_u;
    Quantity {
      { Name Den; Value {Term {[DensityEval[ElementNum[]] ];In CoreAll; Jacobian JacVol; }} }
      { Name u ; Value {
        Term { [Norm[ Vector[{ux},{uy},0]*1e3]]; In CoreAll ; Jacobian JacVol ; } } } // Calculation of deflection in polar system
      { Name ux ; Value {
        Term { [Norm[ Vector[{ux},0,0]*1e3]]; In CoreAll ; Jacobian JacVol ; } } }    // Calculation of deflection in polar system
      { Name uy ; Value {
        Term { [Norm[ Vector[0,{uy},0]*1e3]]; In CoreAll ; Jacobian JacVol ; } } }    // Calculation of deflection in polar system
	  { Name SigMis ; Value {
        Local { [ Norm[Mesies[{d ux},{d uy}] ]]; In CoreAll ; Jacobian JacVol ; } } }	   // Calculation of VonMissies Stress
    }
  }
  /*Defenition of Distance equations*/ 
  { Name distancecalc; NameOfFormulation DistanceCalcuations;
    Quantity {
      { Name Distance; Value {Term {[{dis}];In VolumeAll; Jacobian JacVol; }} }
    }
  }
}

PostOperation {
  { Name Map; NameOfPostProcessing Sys_Mag_2D;
    Operation {
      Print[ b, OnElementsOf CoreAll, File StrCat["b",AngleSTR,".pos"] ];
      Print[ a, OnElementsOf CoreAll, File StrCat["a",AngleSTR,".pos"] ];
      Print[ CurrentStream, OnElementsOf Conductors, File StrCat["stream.pos"] ];
      Print[ AREA, OnElementsOf Conductors, File StrCat["stream.pos"] ];
    }
  }
  { Name Torque; NameOfPostProcessing Sys_Mag_2D;
    Operation {
	  Print[ TorqueMaxwell[AirBoundry], OnGlobal, Format Table, File  StrCat["T",AngleSTR,".pos"],
      SendToServer Sprintf("Output/Magnet %g/Torque 45 degrees [Nm]", 0), Color "Ivory"  ]; 	// Torque by stress tensor
    }
  }
  /* Torqueripple*/ 
  { Name Torque_ripple; NameOfPostProcessing Sys_Mag_2D;
    Operation {
      //Print[ b, OnElementsOf CoreAll, AppendExpressionToFileName $Time,File  StrCat["b_rip_.pos"] ];
	  Print[ TorqueMaxwell[AirBoundry], OnGlobal, AppendExpressionToFileName $Time, Format Table, File  StrCat["T_rip_.pos"],
      SendToServer Sprintf("Output/Magnet %g/Torque 45 degrees [Nm]", 0), Color "Ivory"  ]; 	// Torque by stress tensor
     }
 } 
/*Defecnition of Mechanical Results to return to user*/ 
  { Name pos; NameOfPostProcessing Elast_u_2D;
    Operation {
	  Print[ u, OnElementsOf CoreAll,File  "u.pos"];								// Maximum Deflection Value
	  Print[ Den, OnElementsOf CoreAll,File  "Density.txt"];								// Maximum Deflection Value
	  Print[ SigMis, OnElementsOf CoreAll,File  "SigMis.pos"];  					// Maximum Vonmisies stress
    }
  } 
  /* Distance Terms*/ 
  { Name Dis; NameOfPostProcessing distancecalc;
    Operation {
      //Print[ Distance, OnElementsOf Volume,File  "DistanceReturned.txt"]; // Distances imported
     }
 } 
}