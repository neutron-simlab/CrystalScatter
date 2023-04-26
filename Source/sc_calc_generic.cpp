#include "sc_calc_generic.h"


SC_Calc_GENERIC::SC_Calc_GENERIC()
{
    calc = new SasCalc_GENERIC_calculation();
}


QStringList SC_Calc_GENERIC::guiLayout()
{
    DT( qDebug() << "guiLayout()" );
    // Each element must be in the form "x;y;prompt;type;tooltip;default" with:
    //  x;y     = index in the grid (0,1,2,....)
    //  prompt  = prompting text label left of the inputfield
    //  type    = Selection : "cbs|...|...|..."
    //            Textinput : "txt|len"
    //            Numericals: "inp|frac|min|max|unit"
    //            CheckBox  : "tog"
    //            Infolabel : "lbl"
    //  tooltip = this is the tooltip set to both prompt label and inputfield (optional)
    //  default = the default value (optional)

    // UPDATE März 2023:
    //   Die GUI wird komplett überarbeitet. Alle Informationen werden aus dem ObjName des GUI-Elements
    //   bezogen. Diese Liste hier wird nur noch für die Grenzwerte/Einheiten/Listen genutzt. D.h. die
    //   ersten beiden Felder (x;y) werden nicht mehr genutzt, sollten für eine Übergangszeit aber noch
    //   eindeutig bleiben.
    //   Diese Liste wird in der Consolenversion auch genutzt.

    static QStringList slFCC =
        {
            "N;0;EditBFactor;inp|3|0.01|999;;1",
            "N;0;EditCeff;inp|4|0.01|999;;0.01",
            "N;0;EditQmax;inp|4;Qmax preset from user;2",
            "N;0;EditRadius;inp;Inner radius;20",
            "N;0;EditRadiusi;inp;Outer radius;0",
            "N;0;EditSigma;inp|4;;0.06",
            "N;0;EditDbeta;inp|4|0|360|°;;2",
            "N;0;EditDebyeWaller;inp|3|0|100|nm;;1",
            "N;0;EditPixelNoX;inp;;128",
            "N;0;EditPixelNoY;inp;;128",
            "N;0;EditPixelX;inp|4|0.001|1|m;;0.01",
            "N;0;EditPixelY;inp|4|0.001|1|m;;0.01",
            "N;0;EditDet;inp|4|0.001|100|m;;10",
            "N;0;P1;inp;;0",       // wird in formpq verwendet
            "N;0;Length;inp|4|0|1000;;1",
            "N;0;SigmaL;inp|4;;0.06",
            //"N;0;ShellNo;inp;;0",
            "N;0;reff;inp;;0",
            "N;0;EditCeffcyl;inp|4|0|999;;0",

            "C;0;Ordis;cbs|Gaussian|Exponential|Onsager|Maier-Saupe|Cut-off|Laguerre|z-dir|isotropic|mirrored Gaussian|mirrored Exponential|mirrored Onsager|mirrored Maier-Saupe|mirrored Cut-off|fiber pattern;;isotropic",
            "C;0;ComboBoxInterior;cbs|homogeneous|core + homogeneous sh|core + inhomogeneous sh|multi-shell|myelin;;homogeneous",
            "C;0;ComboBoxParticle;cbs|sphere|cylinder|disk|vesicle|cube|ellipsoid|triaxial ellipsoid|super ellipsoid, barrel|superball|excluded volume chain|Kratky Porod chain;;Sphere",

            "C;0;ComboBoxPeak;cbs|Lorentzian|Gaussian|mod. 1 Lorentzian|mod. 2 Lorentzian|Pseudo-Voigt|Pearson VII|Gamma|Anisotropic Gaussian;;",
            // DEF='Anisotropic Gaussian' sieht blöd aus weil dann 'EditPeakPar' sichtbar ist, aber wegen LType=None diese
            //          ComboBox disabled (=unsichtbar) gesetzt wird.

            "C;0;LType;cbs|Lamellae|hex cyl|sq cyl|rec cyl|BCC|FCC|HCP|SC|BCT|Ia3d|Pn3m|Im3m|None|CP-Layers"
                "|2D-Hex, GISAXS|2D-square, GISAXS|1D-lam, GISAXS|Fd3m, diamond|orthorombic spheres|QC"
                "|Percus-Yevick|Teubner-Strey|Pm3n, A15;;None",

            "N;0;EditPeakPar;inp;;0",
            "N;0;EditAzi;inp;;80",
            //"N;0;RotX;inp|2|0|360|°;;0",
            //"N;0;RotY;inp|2|0|360|°;;0",

            "N;0;EditDomainSize;inp|2|0|1000|nm;Radial domain size in nm;250",
            "N;0;EditRho;inp;;0",
            "N;0;iso;inp|5|0|1000;;1",    // Multiplikator für radintensity (generic)
            "N;0;I0;inp|5|0|100000;;0",
            "N;0;Base;inp|5|-10000|10000;;0",
            "N;0;ifluc;inp|4|0|1000;;0",
            "N;0;rfluc;inp|4|0|1000;;0",
            //---
            "N;0;uca;inp|2|0.01|100|nm;;80",
            "N;0;ucb;inp|2|0.01|100|nm;;21",
            "N;0;ucc;inp|2|0.01|100|nm;;21",
            "N;0;ucalpha;inp|3|0|360|°;;90",
            "N;0;ucbeta;inp|3|0|360|°;;90",
            "N;0;ucgamma;inp|3|0|360|°;;90",
            "N;0;ucpsi;inp|3|0|360|°;;0",
            "N;0;ucn1;inp;;1",
            "N;0;ucn2;inp;;0",
            "N;0;ucn3;inp;;0",
            "N;0;theta;inp|3|0|360|°;;0",  // -> PolTheta   Begriffe werden nicht geändert, damit
            "N;0;phi;inp|3|0|360|°;;0",    // -> PolPhi     die Daten von vorher lesbar bleiben.
            "N;0;rotTheta;inp|3|0|360|°;;0",
            "N;0;rotPhi;inp|3|0|360|°;;0",
            "N;0;EditWavelength;inp|5|0.001|200|nm;;0.154",

            "N;0;acpl;inp;;0", // Editastack
            "N;0;bcpl;inp;;0", // Editbstack

            "N;0;Alfa;inp;;0",

            //---
            "T;0;RadButDebyeScherrer;tog;;0",
            "T;0;RadioButtonPara;tog;;0",
            "T;0;CheckBoxTwinned;tog;;0",
            "T;0;CheckBoxWAXS;tog;;0",

            // Spezielle Einträge, die in die Params-Struktur übernommen werden. Die GUI-Elemente sind aber
            // im Bereich der globalen Eingaben auf der linken Seite zu finden.
            "N;0;SigX;inp|3|0.01|9999;editdom1;40",
            "N;0;SigY;inp|3|0.01|9999;editdom2;40",
            "N;0;SigZ;inp|3|0.01|9999;editdom3;40",

            "N;0;VAx1;inp|3|-1000|1000;Ax1;1",
            "N;0;VAx2;inp|3|-1000|1000;Ax2;0",
            "N;0;VAx3;inp|3|-1000|1000;Ax3;0",
            "N;0;Ay1;inp|3|-1000|1000;Ay1;0",
            "N;0;Ay2;inp|3|-1000|1000;Ay2;1",
            "N;0;Ay3;inp|3|-1000|1000;Ay3;0",
            "N;0;Az1;inp|3|-1000|1000;Az1;0",
            "N;0;Az2;inp|3|-1000|1000;Az2;0",
            "N;0;Az3;inp|3|-1000|1000;Az3;1",

            "N;0;BeamPosX;inp|3|-1000|1000;;0",
            "N;0;BeamPosY;inp|3|-1000|1000;;0",

            "I;0;HKLmax;int|0|1|20;;3",
            "I;0;GridPoints;int|0|16|1024;;64",

            // Bzgl. der Eingabe von QMax gibt es zwei Wertefelder:
            // "G;0;EditQmax;", -> ist oben schon definiert
            "N;0;CalcQmax;inp|4;Qmax calculated from data header above;",
            // Unterschieden werden die beiden Eingaben über die Radiobuttons
            "T;0;EditQmaxData;tog;;",   // auch wenn das Radiobuttons sind
            "T;0;EditQmaxPreset;tog;;", // -"-

        };
    return slFCC;
}

void SC_Calc_GENERIC::prepareData( _dataGetter dg )
{
    DT( qDebug() << "prepareData()" );
    _valueTypes val, val1, val2;

    // Common settings
    (*dg)( "RadioButtonQ1", val );  calc->setRadQ1( val.value > 0 );
    (*dg)( "RadioButtonQ2", val );  calc->setRadQ2( val.value > 0 );
    (*dg)( "RadioButtonQ4", val );  calc->setRadQ4( val.value > 0 );
    (*dg)( "ExpandImage", val );    calc->setExpandImage( val.value > 0 );
    (*dg)( "GridPoints", val );     calc->setGridPoints( val.value );
    (*dg)( "HKLmax", val );         calc->setHKLmax( val.value );

    // Diese A?? sind jetzt einzelne Elemente und werden einzeln von der GUI gelesen!
    // im Main werden diese unter Ax1, Ax2, Ax3 als Vektoren gespeichert.
    // TODO: die UI-Elemente umbenennen... Dann sollten hier wieder die Arrays genutzt werden können.
    (*dg)( "Ax1", val );    calc->setAx1( val.vec );
    (*dg)( "Ax2", val );    calc->setAx2( val.vec );
    (*dg)( "Ax3", val );    calc->setAx3( val.vec );
    (*dg)( "SigXYZ", val );         calc->setSigXYZ( val.vec );

    (*dg)( "BeamPosX", val  );
    (*dg)( "BeamPosY", val1 );      calc->setBeamStop( val.value, val1.value );

    // specific settings
    (*dg)( "EditRadius", val );      calc->setRadiusF( val.value );
    (*dg)( "EditRadiusi", val );     calc->setRadiusI( val.value );
    (*dg)( "EditSigma", val );       calc->setSigmaF( val.value );
    (*dg)( "EditQmax", val );        calc->setQMax( val.value );
    (*dg)( "CheckBoxTwinned", val ); calc->setCheckBoxTwinned( val.checked );
    (*dg)( "Ordis", val );           calc->setOrdis( val.select );
    (*dg)( "EditDebyeWaller", val ); calc->setDisplacement( val.value );
    (*dg)( "EditCeff", val );        calc->setCeffF( val.value );
    (*dg)( "EditBFactor", val );     calc->setBFactorF( val.value );
    (*dg)( "EditDbeta", val );       calc->setDBetaF( val.value );

    (*dg)( "EditCeffcyl", val );     calc->setCeffCyl( val.value );

    (*dg)( "LType", val );            calc->setLType( val.select );

    (*dg)( "ComboBoxInterior", val ); calc->setComboBoxInterior( val.select );
    (*dg)( "ComboBoxParticle", val ); calc->setComboBoxParticle( val.select );
    (*dg)( "ComboBoxPeak", val );     calc->setComboBoxPeak( val.select );
    (*dg)( "EditPeakPar", val );      calc->setPeakPar( val.value );
    (*dg)( "EditAzi", val );          calc->setAzi( val.value );
    (*dg)( "RotX", val  );           //
    (*dg)( "RotY", val1 );           //
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

    (*dg)( "Alfa", val );             calc->setAlphash( val.value );

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

    (*dg)( "theta", val );            calc->setPolTheta( val.value );
    (*dg)( "phi", val );              calc->setPolPhi( val.value );
    (*dg)( "rotTheta", val );         calc->setRotTheta( val.value );
    (*dg)( "rotPhi", val );           calc->setRotPhi( val.value );
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

    //(*dg)( "EditRelDis", val );       calc->setRelDis( val.value );
    //(*dg)( "EditDist", val );         calc->setDist( val.value );
    (*dg)( "EditDomainSize", val );   calc->setDomainsize( val.value );

    (*dg)( "P1", val );               calc->setP1( val.value );
    (*dg)( "SigmaL", val );           calc->setSigmaL( val.value );
    (*dg)( "Length", val );           calc->setLength( val.value );
    //(*dg)( "ShellNo", val );          calc->setShellNo( val.value );
    (*dg)( "reff", val );             calc->setReff( val.value );
    (*dg)( "acpl", val );             calc->setAcpl( val.value );
    (*dg)( "bcpl", val );             calc->setBcpl( val.value );
}
