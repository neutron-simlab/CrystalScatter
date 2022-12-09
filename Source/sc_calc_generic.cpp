#include "sc_calc_generic.h"


SC_Calc_GENERIC::SC_Calc_GENERIC() : SC_Calc()
{
    calc = new SasCalc_GENERIC_calculation();
}


QStringList SC_Calc_GENERIC::guiLayout()
{
    // Each element must be in the form "x;y;prompt;type;tooltip;default" with:  ==> sc_calc.h for details
    //  x;y     = index in the grid (0,1,2,....)
    //  prompt  = prompting text label left of the inputfield
    //  type    = Selection : "cbs|...|...|..."
    //             Fittable : "cbsfit|...|...|..."
    //            Textinput : "txt|len"
    //            Numericals: "inp|frac|min|max|unit"
    //              Fittable: "inpfit|frac|min|max|unit"
    //            CheckBox  : "tog"
    //            Infolabel : "lbl"
    //  tooltip = this is the tooltip set to both prompt label and inputfield (optional)
    //  default = the default value (optional)
    static QStringList slFCC =
    { "0;0;EditBFactor;inpfit|3|0.01|999;;1",
      "0;1;EditCeff;inpfit|4|0.01|999;;0.01",
      "0;2;EditQmax;inpfit|4;;2",
      "0;3;EditRadius;inpfit;Inner radius;20",
      "0;4;EditRadiusi;inpfit;Outer radius;0",
      "0;5;EditSigma;inpfit|4;;0.06",
      "0;6;EditStretch;inp;;1",
      "0;7;EditDbeta;inpfit|4|0|360|°;;2",
      "0;8;EditDebyeWaller;inpfit|3|0|100|nm;;1",
      "X0;9;MaxIter;inp;;1",  // normal :5
      "0;10;EditPixelNoX;inp;;128",
      "0;11;EditPixelNoY;inp;;128",
      "0;12;EditPixelX;inp|4|0.001|1|m;;0.01",
      "0;13;EditPixelY;inp|4|0.001|1|m;;0.01",
      "0;14;EditDet;inp|4|0.001|100|m;;10",
      "0;15;P1;inp;;0",
      "0;16;Length;inpfit|4|0|1000;;1",
      "0;17;SigmaL;inpfit|4;;0.06",
      "0;18;ShellNo;inp;;0",
      "0;19;reff;inp;;0",

      //---
      "1;0;Ordis;cbs|Gaussian|Exponential|Onsager|Maier-Saupe|Cut-off|Laguerre|z-dir|isotropic|mirrored Gaussian|mirrored Exponential|mirrored Onsager|mirrored Maier-Saupe|mirrored Cut-off|fiber pattern|TEST:(x²*y³);;isotropic",
      "1;1;ComboBoxInterior;cbs|homogeneous|core + homogeneous sh|core + inhomogeneous sh;;homogeneous",
      "1;2;ComboBoxParticle;cbs|Sphere|Cylinder|Disk|Vesicle|Cube|Ellipsoid|Triaxial ellipsoid|Super ellipsoid;;Sphere",
      "1;3;ComboBoxPeak;cbs|Lorentzian|Gaussian|mod. 1 Lorentzian|mod. 2 Lorentzian|Pseudo-Voigt|Pearson VII|Gamma|Anisotropic Gaussian;;Anisotropic Gaussian",
      "1;4;EditPeakPar;inp;;0",
      "1;5;EditAzi;inpfit;;80",
      "1;6;Rot_X;inp|2|0|360|°;;0",
      "1;7;Rot_Y;inp|2|0|360|°;;0",
      "X1;8;Rot_Z;inp|2|0|360|°;;0",
      "X1;9;Rot_Angle;inp|2|0|360|°;;0",
      "X1;10;Tilt_Angle;inp|2|0|360|°;;0",
      "X1;11;Editx;inp;;0",
      "X1;12;Edity;inp;;0",
      "X1;13;EditAnglexy;inp;;0",
            // Achtung: wenn oben das X wegfällt, müssen die beiden nächsten Werte verschoben werden
      "1;12;EditRelDis;inp|2;;1",
      "1;13;EditDist;inp|2;;1",
      "1;14;EditDomainSize;inpfit|2|0|1000|nm;Radial domain size in nm;250",
      "1;15;EditRho;inp;;0",
      "1;16;iso;inpfit|5|0|1000;;1",    // Multiplikator für radintensity (generic)
      "1;17;I0;inpfit|5|0.01|999999;;1000",
      "1;18;Base;inpfit|5|-10000|10000;;0",
      "1;19;ifluc;inp|4|0|1000;;0",
      "1;20;rfluc;inp|4|0|1000;;0",
      //---
      "2;0;unit cell definiton;lbl;;",
      "2;1;uca;inpfit|2|0.01|100|nm;;80",
      "2;2;ucb;inpfit|2|0.01|100|nm;;21",
      "2;3;ucc;inpfit|2|0.01|100|nm;;21",
      "2;4;ucalpha;inp|3|0|360|°;;90",
      "2;5;ucbeta;inp|3|0|360|°;;90",
      "2;6;ucgamma;inp|3|0|360|°;;90",
      "2;7;ucpsi;inpfit|3|0|360|°;;0",
      "2;8;ucn1;inpfit;;1",
      "2;9;ucn2;inpfit;;0",
      "2;10;ucn3;inpfit;;0",
      "2;11;theta;inpfit|3|0|360|°;;0",
      "2;12;phi;inpfit|3|0|360|°;;0",
      "2;13;EditWavelength;inp|5|0.001|200|nm;;0.154",
      "X2;14;EditWAXSangle;inp;;90",

      "2;15;acpl;inp;;0",
      "2;16;bcpl;inp;;0",
      "2;17;por;inp;;0",

      "2;18;LType;cbs|Lam|hex cyl|sq cyl|rec cyl|BCC|FCC|HCP|SC|BCT|Ia3d|Pn3m|Im3m|None|CP-Layers"
            "|2D-Hex, GISAXS|2D-square, GISAXS|1D-lam, GISAXS|Fd3m, diamond|orthorombic spheres|QC;;None",
      //if ComboBoxLattice.ItemIndex=0 then ltype:=0;   (* Lam *)
      //if ComboBoxLattice.ItemIndex=1 then ltype:=1;   (* hex cyl *)
      //if ComboBoxLattice.ItemIndex=2 then ltype:=2;   (* sq Cyl *)
      //if ComboBoxLattice.ItemIndex=3 then ltype:=3;   (* rec cyl *)
      //if ComboBoxLattice.ItemIndex=4 then ltype:=4;   (* BCC *)
      //if ComboBoxLattice.ItemIndex=5 then ltype:=5;   (* FCC *)
      //if ComboBoxLattice.ItemIndex=6 then ltype:=6;   (* HCP *)
      //if ComboBoxLattice.ItemIndex=7 then ltype:=7;   (* SC *)
      //if ComboBoxLattice.ItemIndex=8 then ltype:=8;   (* BCT *)
      //if ComboBoxLattice.ItemIndex=9 then ltype:=9;   (* Ia3d *)
      //if ComboBoxLattice.ItemIndex=10 then ltype:=10; (* Pn3m *)
      //if ComboBoxLattice.ItemIndex=11 then ltype:=11; (* Im3m *)
      //if ComboBoxLattice.ItemIndex=12 then ltype:=12;  (* none *)
      //if ComboBoxLattice.ItemIndex=13 then ltype:=13;  (* CP-layers *)
      //if ComboBoxLattice.ItemIndex=14 then ltype:=14;  (* 2D-Hex, GISAXS *)
      //if ComboBoxLattice.ItemIndex=15 then ltype:=15;  (* 2D-square, GISAXS *)
      //if ComboBoxLattice.ItemIndex=16 then ltype:=16;  (* 1D-lam, GISAXS *)
      //if ComboBoxLattice.ItemIndex=17 then ltype:=17;  (* Fd3m, diamond *)
      //if ComboBoxLattice.ItemIndex=18 then ltype:=18;  (* orthorombic spheres *)
      //if ComboBoxLattice.ItemIndex=19 then ltype:=19;  (* QC *)

      //---
      "X3;0;CheckBoxf2q;tog;;0",
      "X3;1;RadioButtonCHS;tog;;0",
      "X3;2;RadioButtonCS;tog;;0",
      "3;3;RadButDebyeScherrer;tog;;0",
      "3;4;RadioButtonPara;tog;;0",
      "X3;5;RadioButtonSolid;tog;;0",
      //"3;6;RadioButtonVertical;tog;;0",
      "3;7;CheckBoxTwinned;tog;;0",
      "3;9;CheckBoxWAXS;tog;;0",

//      "3;11;CheckBox10;tog;;0",

      // Spezielle Einträge, die in die Params-Struktur übernommen werden. Die GUI-Elemente sind aber
      // im Bereich der globalen Eingaben auf der linken Seite zu finden.
      "G;0;Editdom1;inpfit|3|0.01|9999;;40",
      "G;0;Editdom2;inpfit|3|0.01|9999;;40",
      "G;0;Editdom3;inpfit|3|0.01|9999;;40",

      // Spezielle Codierung für gesperrte Elemente aus dem generischen Teil
      "X;Uvec;Vvec;Nvec",

        /*
    params.amax = 10; // TODO, noch kein GUI-Element vorhanden
    params.bmax = 10;   --> mehr für die Grafik in Scatter
    params.cmax = 10;
         */

    };
    return slFCC;
}

void SC_Calc_GENERIC::prepareData( _dataGetter dg )
{
    _valueTypes val, val1, val2;

    // Common settings
    (*dg)( "RadioButtonQ1", val );  calc->setRadQ1( val.value > 0 );
    (*dg)( "RadioButtonQ2", val );  calc->setRadQ2( val.value > 0 );
    (*dg)( "RadioButtonQ4", val );  calc->setRadQ4( val.value > 0 );
    (*dg)( "ExpandImage", val );    calc->setExpandImage( val.value > 0 );
    (*dg)( "EditGridPoints", val ); calc->setGridPoints( val.value );
    (*dg)( "Edithklmax", val );     calc->setHKLmax( val.value );
    //xx//(*dg)( "Uvec", val );           calc->setUvec( val.vec );
    //xx//(*dg)( "Vvec", val );           calc->setVvec( val.vec );
    //xx//(*dg)( "Nvec", val );           calc->setNvec( val.vec );
    (*dg)( "Ax1", val );            calc->setAx1( val.vec );
    (*dg)( "Ax2", val );            calc->setAx2( val.vec );
    (*dg)( "Ax3", val );            calc->setAx3( val.vec );
    (*dg)( "SigXYZ", val );         calc->setSigXYZ( val.vec );

    (*dg)( "BeamPosX", val  );
    (*dg)( "BeamPosY", val1 );      calc->setBeamStop( val.value, val1.value );

    // specific settings
    (*dg)( "EditRadius", val );      calc->setRadiusF( val.value );
    (*dg)( "EditSigma", val );       calc->setSigmaF( val.value );
    (*dg)( "EditQmax", val );        calc->setQMax( val.value );
    (*dg)( "EditStretch", val );     calc->setStretch( val.value );
    (*dg)( "CheckBoxTwinned", val ); calc->setCheckBoxTwinned( val.checked );
    (*dg)( "Ordis", val );           calc->setOrdis( val.select );
    (*dg)( "EditDebyeWaller", val ); calc->setDisplacement( val.value );
    //xx//(*dg)( "MaxIter", val );         calc->setMaxIter( val.value );
    (*dg)( "EditCeff", val );        calc->setCeffF( val.value );
    (*dg)( "EditBFactor", val );     calc->setBFactorF( val.value );
    (*dg)( "EditDbeta", val );       calc->setDBetaF( val.value );

    (*dg)( "LType", val );            calc->setLType( val.select );

    (*dg)( "ComboBoxInterior", val ); calc->setComboBoxInterior( val.select );
    (*dg)( "ComboBoxParticle", val ); calc->setComboBoxParticle( val.select );
    (*dg)( "ComboBoxPeak", val );     calc->setComboBoxPeak( val.select );
    (*dg)( "EditPeakPar", val );      calc->setPeakPar( val.value );
    (*dg)( "EditAzi", val );          calc->setAzi( val.value );
    (*dg)( "Rot_X", val  );           //
    (*dg)( "Rot_Y", val1 );           //
    //xx//(*dg)( "Rot_Z", val2 );           calc->setRotation( Double3(val.value,
    //xx//                                                             val1.value,
    //xx//                                                             val2.value) );
    //xx//(*dg)( "Rot_Angle", val );        calc->setRotAngle( val.value );
    //xx//(*dg)( "Tilt_Angle", val );       calc->setTiltAng( val.value );
    //xx//(*dg)( "Editx", val  );           //
    //xx//(*dg)( "Edity", val1 );           calc->setxycur( val.value,
    //xx//                                                  val1.value );
    //xx//(*dg)( "EditAnglexy", val );      calc->setanglecur( val.value );
    (*dg)( "CheckBoxWAXS", val );     calc->setCheckBoxWAXS( val.checked );
    //xx//(*dg)( "CheckBoxf2q", val );      calc->setCheckBoxf2q( val.checked );
    //(*dg)( "CheckBoxLeftCircle", val );  calc->setCheckBoxLeftCircle( val.checked );
    //(*dg)( "CheckBoxLeftHoriz", val );   calc->setCheckBoxLeftHoriz( val.checked );
    //(*dg)( "CheckBoxLeftVert", val );    calc->setCheckBoxLeftVert( val.checked );
    //xx//(*dg)( "RadioButtonCHS", val );   calc->setRadioButtonCHS( val.checked );
    //xx//(*dg)( "RadioButtonCS", val );    calc->setRadioButtonCS( val.checked );
    (*dg)( "RadButDebyeScherrer", val ); calc->setRadioButtonDebyeScherrer( val.checked );
    (*dg)( "RadioButtonPara", val );  calc->setRadioButtonPara( val.checked );
    //xx//(*dg)( "RadioButtonSolid", val ); calc->setRadioButtonSolid( val.checked );
    //(*dg)( "RadioButtonVertical", val ); calc->setRadioButtonVertical( val.checked );

//    (*dg)( "CheckBox10", val );     calc->setCheckBox10( val.checked );

    //{NV} - unit cell definiton
    (*dg)( "uca", val );              calc->setUCA( val.value );
    (*dg)( "ucb", val );              calc->setUCB( val.value );
    (*dg)( "ucc", val );              calc->setUCC( val.value );
    (*dg)( "ucalpha", val );          calc->setUCalpha( val.value );
    (*dg)( "ucbeta", val );           calc->setUCbeta( val.value );
    (*dg)( "ucgamma", val );          calc->setUCgamma( val.value );
    (*dg)( "ucpsi", val );            calc->setUCpsi( val.value );
    (*dg)( "ucn1", val );             calc->setUCn1( val.value );
    (*dg)( "ucn2", val );             calc->setUCn2( val.value );
    (*dg)( "ucn3", val );             calc->setUCn3( val.value );

    (*dg)( "theta", val );            calc->setTheta( val.value );
    (*dg)( "phi", val );              calc->setPhi( val.value );
    (*dg)( "EditRho", val );          calc->setRho( val.value );

    (*dg)( "iso", val );              calc->setIso( val.value );
    (*dg)( "I0", val );               calc->setIZero( val.value );
    (*dg)( "Base", val );             calc->setBase( val.value );
    (*dg)( "ifluc", val );            calc->setIFluc( val.value );
    (*dg)( "rfluc", val );            calc->setRFluc( val.value );

    (*dg)( "EditPixelNoX", val );     calc->setpixnox( val.value );
    (*dg)( "EditPixelNoY", val );     calc->setpixnoy( val.value );
    (*dg)( "EditPixelX", val );       calc->setpixx( val.value );
    (*dg)( "EditPixelY", val );       calc->setpixy( val.value );
    (*dg)( "EditDet", val );          calc->setdet( val.value );
    //xx//(*dg)( "EditWAXSangle", val );    calc->setxrdalf( val.value );
    (*dg)( "EditWavelength", val );   calc->setwave( val.value );

    (*dg)( "EditRelDis", val );       calc->setRelDis( val.value );
    (*dg)( "EditDist", val );         calc->setDist( val.value );
    (*dg)( "EditDomainSize", val );   calc->setDomainsize( val.value );

    (*dg)( "P1", val );               calc->setP1( val.value );
    (*dg)( "SigmaL", val );           calc->setSigmaL( val.value );
    (*dg)( "Length", val );           calc->setLength( val.value );
    (*dg)( "ShellNo", val );          calc->setShellNo( val.value );
    (*dg)( "reff", val );             calc->setReff( val.value );
    (*dg)( "por", val );              calc->setPor( val.value );
    (*dg)( "acpl", val );             calc->setAcpl( val.value );
    (*dg)( "bcpl", val );             calc->setBcpl( val.value );

    // Jetzt noch die Lattice-Parameter aus der GUI übertragen (speziell für den 2d-Fit)
    (*dg)( "LATTcols", val );   int cols        = val.value;
    (*dg)( "LATTrows", val );   int rows        = val.value;
    (*dg)( "LATTcenx", val );   double centerx  = val.value;
    (*dg)( "LATTceny", val );   double centery  = val.value;
    (*dg)( "LATTwlen", val );   double wavelen  = val.value;
    (*dg)( "LATTdist", val );   double distance = val.value;
    (*dg)( "LATTpixx", val );   double pixx     = val.value / 1000.0;   // Anzeige in mm
    (*dg)( "LATTpixy", val );   double pixy     = val.value / 1000.0;   // Nutzung in m (da so im normalen Calc auch)
    calc->setLattPar( cols, rows, centerx, centery, wavelen, distance, pixx, pixy );

}
