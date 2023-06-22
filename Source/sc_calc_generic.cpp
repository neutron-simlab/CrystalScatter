#include "sc_calc_generic.h"


SC_Calc_GENERIC::SC_Calc_GENERIC()
{
    calc = new SasCalc_GENERIC_calculation();
}


QStringList SC_Calc_GENERIC::guiLayoutNeu()
{
    DT( qDebug() << "guiLayout()" );

    // UPDATE März 2023:
    //   Die GUI wird komplett überarbeitet. Alle Informationen werden aus dem ObjName des GUI-Elements
    //   bezogen. Diese Liste hier wird nur noch für die Grenzwerte/Einheiten/Listen genutzt. Im ersten
    //   Feld (x) steht jetzt auch eine Datentypskennung, das zweite Feld (y) wird nicht mehr genutzt.
    //   Diese Liste wird in der Consolenversion auch genutzt.

    // Each element must be in the form "type;kenn;typespec;tooltip;default" with:
    //  type     = C:Selection, N:Double, I:Integer, T:Toggle, O:DoubleOut
    //  kenn     = internal parameter name to connect to correct gui element
    //  typespec = C:Selection : "...|...|..."  (required)
    //             N:Double    : "frac|min|max|unit"  (optional, default: "2|-10000|+10000|")
    //             I:Integer   : "min|max|unit"  (optional, default: "-10000|+10000|")
    //             T:Toggle    : (empty)
    //             O:DoubleOut : (empty)
    //  tooltip  = this is the tooltip set to both prompt label and inputfield (optional)
    //  default  = the default value (optional)

    static QStringList slGUI =
        {
            // Lattice group
            "C;LType;Lamellae|hex cyl|sq cyl|rec cyl|BCC|FCC|HCP|SC|BCT|Ia3d|Pn3m|Im3m|None|CP-Layers"
                "|2D-Hex, GISAXS|2D-square, GISAXS|1D-lam, GISAXS|Fd3m, diamond|orthorombic spheres|QC"
                "|Percus-Yevick|Teubner-Strey|Pm3n, A15;;None",
            "N;uca;2|0.01|100|nm;;80",
            "N;ucb;2|0.01|100|nm;;21",
            "N;ucc;2|0.01|100|nm;;21",
            "N;ucalpha;3|0|360|°;;90",
            "N;ucbeta;3|0|360|°;;90",
            "N;ucgamma;3|0|360|°;;90",
            "N;EditCeff;4|0.01|999;;0.01",
            "T;CheckBoxTwinned;;;0",
            "N;EditCeffcyl;4|0|999;;0",
            "N;reff;;;0",
            "N;acpl;;;0", // Editastack
            "N;bcpl;;;0", // Editbstack
            "N;ifluc;4|0|1000;;0",
            "N;rfluc;4|0|1000;;0",
            "O;EditRelDis;;;",
            "O;EditDist;;;",

            // Particle group
            "C;ComboBoxParticle;sphere|cylinder|disk|vesicle|cube|ellipsoid|triaxial ellipsoid|super ellipsoid, barrel"
                "|superball|excluded volume chain|Kratky Porod chain;;Sphere",
            "C;Ordis;Gaussian|Exponential|Onsager|Maier-Saupe|Cut-off|Laguerre|z-dir|isotropic|mirrored Gaussian"
                "|mirrored Exponential|mirrored Onsager|mirrored Maier-Saupe|mirrored Cut-off|fiber pattern;;isotropic",
            "C;ComboBoxInterior;homogeneous|core + homogeneous sh|core + inhomogeneous sh|multi-shell|myelin;;homogeneous",
            "N;EditRadius;;Inner radius;20",
            "N;EditRadiusi;;Outer radius;0",
            "N;EditSigma;4;;0.06",
            "N;EditDbeta;4|0|360|°;;2",
            "N;Length;4|0|1000;;1",
            "N;SigmaL;4;;0.06",
            "N;Alpha;;Internal 'alphash';0",
            "N;EditRho;;;0",
            "N;RotAlpha;;Internal 'alpha';0",

            // Peak Shape group
            "C;ComboBoxPeak;Lorentzian|Gaussian|mod. 1 Lorentzian|mod. 2 Lorentzian|Pseudo-Voigt|Pearson VII|Gamma|Anisotropic Gaussian;;",
            "N;EditPeakPar;;;0",
            "N;EditDebyeWaller;3|0|100|nm;;1",
            "N;EditAzi;;;80",
            "N;EditDomainSize;2|0|1000|nm;Radial domain size in nm;250",
            "T;RadButDebyeScherrer;;;0",
            "T;RadioButtonPara;;;0",
            "N;VAx1;3|-1000|1000;Ax1;1",
            "N;VAx2;3|-1000|1000;Ax2;0",
            "N;VAx3;3|-1000|1000;Ax3;0",
            "N;Ay1;3|-1000|1000;Ay1;0",
            "N;Ay2;3|-1000|1000;Ay2;1",
            "N;Ay3;3|-1000|1000;Ay3;0",
            "N;Az1;3|-1000|1000;Az1;0",
            "N;Az2;3|-1000|1000;Az2;0",
            "N;Az3;3|-1000|1000;Az3;1",
            "N;SigX;3|0.01|9999;editdom1;40",
            "N;SigY;3|0.01|9999;editdom2;40",
            "N;SigZ;3|0.01|9999;editdom3;40",

            // Orientation group
            "N;ucpsi;3|0|360|°;;0",
            "N;ucn1;;;1",
            "N;ucn2;;;0",
            "N;ucn3;;;0",
            "N;theta;3|0|360|°;;0",  // -> PolTheta   Begriffe werden nicht geändert, damit
            "N;phi;3|0|360|°;;0",    // -> PolPhi     die Daten von vorher lesbar bleiben.
            "N;rotTheta;3|0|360|°;;0",
            "N;rotPhi;3|0|360|°;;0",

            // Calculation group
            "I;GridPoints;16|2049;;64",
            "I;HKLmax;1|20;;3",
            "N;EditQmax;4;Qmax preset from user;2",                   // Das genutzte QMax wird über die
            "N;CalcQmax;4;Qmax calculated from data header above;",   // Radiobuttons darunter ausgewählt
            "T;EditQmaxData;;;",   // auch wenn das Radiobuttons sind
            "T;EditQmaxPreset;;;", // -"-

            // Experiment group
            "N;EditPixelNoX;;;128",
            "N;EditPixelNoY;;;128",
            "N;EditPixelX;4|0.001|1|m;;0.01",
            "N;EditPixelY;4|0.001|1|m;;0.01",
            "N;EditDet;4|0.001|100|m;;10",
            "N;EditWavelength;5|0.001|200|nm;;0.154",
            "N;BeamPosX;3|-1000|1000;;0",
            "N;BeamPosY;3|-1000|1000;;0",

            // Controls group
            "N;EditBFactor;3|0.01|999;;1",
            "T;CheckBoxWAXS;;;0",
            "N;P1;;;0",       // wird in formpq verwendet

            // Pixel Manipulation group
            "N;iso;5|0|1000;;1",    // Multiplikator für radintensity (generic)
            "N;I0;5|0|100000;;0",
            "N;Base;5|-10000|10000;;0",
        };
    return slGUI;
}

void SC_Calc_GENERIC::prepareData( _dataGetter dg )
{
    DT( qDebug() << "prepareData()" );
    _valueTypes val, val1;

    // Lattice group
    (*dg)( "LType", val );               calc->setLType( val.select );
    (*dg)( "uca", val );                 calc->setUCA( val.value );
    (*dg)( "ucb", val );                 calc->setUCB( val.value );
    (*dg)( "ucc", val );                 calc->setUCC( val.value );
    (*dg)( "ucalpha", val );             calc->setUCalpha( val.value );
    (*dg)( "ucbeta", val );              calc->setUCbeta( val.value );
    (*dg)( "ucgamma", val );             calc->setUCgamma( val.value );
    (*dg)( "EditCeff", val );            calc->setCeffF( val.value );
    (*dg)( "CheckBoxTwinned", val );     calc->setCheckBoxTwinned( val.checked );
    (*dg)( "EditCeffcyl", val );         calc->setCeffCyl( val.value );
    (*dg)( "reff", val );                calc->setReff( val.value );
    (*dg)( "acpl", val );                calc->setAcpl( val.value );
    (*dg)( "bcpl", val );                calc->setBcpl( val.value );
    (*dg)( "ifluc", val );               calc->setIFluc( val.value );
    (*dg)( "rfluc", val );               calc->setRFluc( val.value );

    // Particle group
    (*dg)( "ComboBoxParticle", val );    calc->setComboBoxParticle( val.select );
    (*dg)( "Ordis", val );               calc->setOrdis( val.select );
    (*dg)( "ComboBoxInterior", val );    calc->setComboBoxInterior( val.select );
    (*dg)( "EditRadius", val );          calc->setRadiusF( val.value );
    (*dg)( "EditRadiusi", val );         calc->setRadiusI( val.value );
    (*dg)( "EditSigma", val );           calc->setSigmaF( val.value );
    (*dg)( "EditDbeta", val );           calc->setDBetaF( val.value );
    (*dg)( "Length", val );              calc->setLength( val.value );
    (*dg)( "SigmaL", val );              calc->setSigmaL( val.value );
    (*dg)( "Alpha", val );               calc->setAlpha( val.value );
    (*dg)( "EditRho", val );             calc->setRho( val.value );
    (*dg)( "RotAlpha", val );            calc->setRotAlpha( val.value );

    // Peak Shape group
    (*dg)( "ComboBoxPeak", val );        calc->setComboBoxPeak( val.select );
    (*dg)( "EditPeakPar", val );         calc->setPeakPar( val.value );
    (*dg)( "EditDebyeWaller", val );     calc->setDisplacement( val.value );
    (*dg)( "EditAzi", val );             calc->setAzi( val.value );
    (*dg)( "EditDomainSize", val );      calc->setDomainsize( val.value );
    (*dg)( "RadButDebyeScherrer", val ); calc->setRadioButtonDebyeScherrer( val.checked );
    (*dg)( "RadioButtonPara", val );     calc->setRadioButtonPara( val.checked );
    (*dg)( "Ax1", val );                 calc->setAx1( val.vec );
    (*dg)( "Ax2", val );                 calc->setAx2( val.vec );
    (*dg)( "Ax3", val );                 calc->setAx3( val.vec );
    (*dg)( "SigXYZ", val );              calc->setSigXYZ( val.vec );

    // Orientation group
    (*dg)( "ucpsi", val );               calc->setUCpsi( val.value );
    (*dg)( "ucn1", val );                calc->setUCn1( val.value );
    (*dg)( "ucn2", val );                calc->setUCn2( val.value );
    (*dg)( "ucn3", val );                calc->setUCn3( val.value );
    (*dg)( "theta", val );               calc->setPolTheta( val.value );
    (*dg)( "phi", val );                 calc->setPolPhi( val.value );
    (*dg)( "rotTheta", val );            calc->setRotTheta( val.value );
    (*dg)( "rotPhi", val );              calc->setRotPhi( val.value );

    // Calculation group
    (*dg)( "GridPoints", val );          calc->setGridPoints( val.value );
    (*dg)( "HKLmax", val );              calc->setHKLmax( val.value );
    (*dg)( "EditQmax", val );            calc->setQMax( val.value );
    (*dg)( "RadioButtonQ1", val );       calc->setRadQ1( val.value > 0 );
    (*dg)( "RadioButtonQ2", val );       calc->setRadQ2( val.value > 0 );
    (*dg)( "RadioButtonQ4", val );       calc->setRadQ4( val.value > 0 );
    (*dg)( "ExpandImage", val );         calc->setExpandImage( val.value > 0 );

    // Experiment group
    (*dg)( "EditPixelNoX", val );        calc->setpixnox( val.value );
    (*dg)( "EditPixelNoY", val );        calc->setpixnoy( val.value );
    (*dg)( "EditPixelX", val );          calc->setpixx( val.value );
    (*dg)( "EditPixelY", val );          calc->setpixy( val.value );
    (*dg)( "EditDet", val );             calc->setdet( val.value );
    (*dg)( "EditWavelength", val );      calc->setwave( val.value );
    (*dg)( "BeamPosX", val  );
    (*dg)( "BeamPosY", val1 );           calc->setBeamStop( val.value, val1.value );

    // Controls group
    (*dg)( "EditBFactor", val );         calc->setBFactorF( val.value );
    (*dg)( "CheckBoxWAXS", val );        calc->setCheckBoxWAXS( val.checked );
    (*dg)( "P1", val );                  calc->setP1( val.value );

    // Pixel Manipulation group
    (*dg)( "iso", val );                 calc->setIso( val.value );
    (*dg)( "I0", val );                  calc->setIZero( val.value );
    (*dg)( "Base", val );                calc->setBase( val.value );
}


void SC_Calc_GENERIC::updateOutputData( _dataSetter ds )
{
    DT( qDebug() << "+++ updateOutputData() +++" );
    _valueTypes val;

    val.value = calc->getRelDis(); (*ds)( "EditRelDis", val );
    val.value = calc->getDist();   (*ds)( "EditDist", val );
}
