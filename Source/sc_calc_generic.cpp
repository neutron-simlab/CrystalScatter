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
    //             Zweites Zeichen ist 'F' für Fittable oder '-' für nicht fittable
    //  kenn     = internal parameter name to connect to correct gui element
    //  typespec = C:Selection : "...|...|..."  (required)
    //             N:Double    : "frac|min|max|unit"  (optional, default: "2|-10000|+10000|")
    //             I:Integer   : "min|max|unit"  (optional, default: "-10000|+10000|")
    //             T:Toggle    : (empty)
    //             O:DoubleOut : (empty)
    //  tooltip  = this is the tooltip set to both prompt label and inputfield (optional)
    //  default  = the default value (optional)

    // Default-Werte vom C:\SimLab\sas-crystal\Scatter Jun 2023\crystal3d1.pas
    // Mit Ausnahme von LType, der bleibt auf None.
    // Nicht direkt gefundene Werte werden mit "//??" markiert.

    static QStringList slGUI =
        {
            // Lattice group, 07.07.23: Bezeichnungen angepasst und zwei neue Werte dazu!
            "C-;LType;Lamellae|Hexagonally packed cylinders (P6/mm)|Square packed cylinders (P4/mm)|Rectangular centered cylinders (cmm)"
                "|BCC (lm3m)|FCC (Fm3m)|HCP (P6/mmc)|SC (Pm3m)|BCT (l4/mm)|Gyroid (Ia3d)|OBDD (Pn3m)|Plumbers Nightmare (Im3m)|None|CP-Layers"
                "|2DHex (GISAXS)|2DSquare (GISAXS)|1DLam (GISAXS)|Fd3m|Orthorombic|LDQ12;Lattice type selection;None",
                //20240301 - weg: "|Percus-Yevick|Teubner-Strey|Pm3n (A15)|P42/mnm (sigma)|Fddd-network"
            "NF;uca;2|0.01|1000| nm;Unit cell dimension a;30",
            "NF;ucb;2|0.01|1000| nm;Unit cell dimension b;30",
            "NF;ucc;2|0.01|1000| nm;Unit cell dimension c;30",
            "NF;ucalpha;3|0|360|°;Unit cell rotation alpha;90",  // Default auf 90°
            "NF;ucbeta;3|0|360|°;Unit cell rotation beta;90",    //
            "NF;ucgamma;3|0|360|°;Unit cell rotation gamma;90",  //
            "NF;EditCeff;4|0|1;Effective volume fraction;1",
            "T-;CheckBoxTwinned;;;0",  //??
            "NF;EditCeffcyl;4|0|999;;0",  //??
            "NF;reff;;Effective radius;0",
            "NF;acpl;5|0|1;Stacking probability A(s);0",
            "NF;bcpl;5|0|1;Stacking probability B(s);0",
            "NF;ifluc;4|0|1000;;0",  //??
            "NF;rfluc;4|0|1000;;0",  //??
            "O-;EditRelDis;;;",
            "O-;EditDist;;;",

            // Particle group
            "C-;ComboBoxParticle;sphere|cylinder|disk|vesicle|cube|ellipsoid|triaxial ellipsoid|super ellipsoid, barrel"
                "|superball|excluded volume chain|Kratky Porod chain;Particle shape selection;Sphere",
            "C-;Ordis;Gaussian|Exponential|Onsager|Maier-Saupe|Cut-off|Laguerre|z-dir|isotropic|mirrored Gaussian"
                "|mirrored Exponential|mirrored Onsager|mirrored Maier-Saupe|mirrored Cut-off|fiber pattern;Orientational distribution functions;isotropic",
            "C-;ComboBoxInterior;homogeneous|core + homogeneous sh|core + inhomogeneous sh|multi-shell|myelin;Particle nature;homogeneous",
            "NF;EditRadius;4|0|1000|nm;Inner radius;1.0",
            "NF;EditRadiusi;4|0|1000|nm;Outer radius;0.5",  //??
            "NF;EditSigma;4;Relative standard deviation;0.07",
            "NF;EditDbeta;4|0|360|°;Average angle;0.4",
            "NF;Length;4|0|1000;Structure length;1",
            "NF;SigmaL;4;Relative standard deviation of the length;0.06",
            "NF;Alpha;;Internal 'alphash';0",  //??
            "NF;EditRho;;;0",  //??
            "NF;RotAlpha;;Internal 'alpha';0",  //??

            // Peak Shape group
            "C-;ComboBoxPeak;Lorentzian|Gaussian|mod. 1 Lorentzian|mod. 2 Lorentzian|Pseudo-Voigt|Pearson VII|Gamma|Anisotropic Gaussian;Select the peak shape;Gaussian",
            "NF;EditPeakPar;4|0|1e6;Peak shape parameter;0",
            "NF;EditDebyeWaller;3|0|100| nm;Average displacement from ideal lattice point;2",
            "NF;EditAzi;;;180",
            "NF;EditDomainSize;2|0|5000| nm;Radial domain size;400",
            "T-;RadButDebyeScherrer;;;0",  //??
            "T-;RadioButtonPara;;;0",  //??
            "N-;VAx1;3|-1000|1000;Ax1;1",  //??
            "N-;VAx2;3|-1000|1000;Ax2;0",
            "N-;VAx3;3|-1000|1000;Ax3;0",
            "N-;Ay1;3|-1000|1000;Ay1;0",
            "N-;Ay2;3|-1000|1000;Ay2;1",
            "N-;Ay3;3|-1000|1000;Ay3;0",
            "N-;Az1;3|-1000|1000;Az1;0",
            "N-;Az2;3|-1000|1000;Az2;0",
            "N-;Az3;3|-1000|1000;Az3;1",
            "NF;SigX;3|0.01|9999;Average domain size X;40",  //?? editdom1,2,3
            "NF;SigY;3|0.01|9999;Average domain size Y;40",
            "NF;SigZ;3|0.01|9999;Average domain size Z;40",

            // Orientation group
            "NF;ucpsi;3|0|360|°;;0",
            "NF;ucn1;;;1",
            "NF;ucn2;;;1",
            "NF;ucn3;;;1",
            "NF;theta;3|0|360|°;;90",  // -> PolTheta   Begriffe werden nicht geändert, damit
            "NF;phi;3|0|360|°;;0",     // -> PolPhi     die Daten von vorher lesbar bleiben.
            "NF;rotTheta;3|0|360|°;;90",
            "NF;rotPhi;3|0|360|°;;90",

            // Calculation group
            "I-;GridPoints;16|2049;Half of the size of each image dimension;64",  //?? - werden nicht auf 20 runtergesetzt
            "I-;HKLmax;1|20;Number of iterations in the h,k,l-loops (Miller indices);5",
            "N-;EditQmax;4|0.001|100| nm-1;Qmax preset from user;1",                   // Das genutzte QMax wird über die
            "N-;CalcQmax;4|0.001|100| nm-1;Qmax calculated from data header above;",   // Radiobuttons darunter ausgewählt
            "T-;EditQmaxData;;Use the Qmax from the data;",   // auch wenn das Radiobuttons sind
            "T-;EditQmaxPreset;;Use the Qmax provided here;", // -"-
            "N-;EditQmin;4|0.0|100| nm-1;Qmin default for 1d graphs;0",
            "I-;EditQsteps;10|10000;Number of steps from Qmin to Qmax for 1d graphs;200",

            // Experiment group
            "NF;EditPixelNoX;1|16|10000;Number of horizontal detector pixel;2048",
            "NF;EditPixelNoY;1|16|10000;Number of vertical detector pixel;2048",
            "NF;EditPixelX;5|0.1|50| mm;Width of one detector pixel;1",
            "NF;EditPixelY;5|0.1|50| mm;Height of one detector pixel;1",
            "NF;EditDet;4|0.001|100| m;Distance Sample - Detector;1",
            "NF;EditWavelength;5|0.001|200| nm;Wavelength;0.09499",
            "NF;BeamPosX;3|-1000|1000;Horizontal center of the beam in pixel coordinates;0",
            "NF;BeamPosY;3|-1000|1000;Vertical center of the beam in pixel coordinates;0",
            "T-;CenterBeam;;Use beam position;",          // Radiobutton für die Beamposition
            "T-;CenterMidpoint;;Use center point;",       // Radiobutton für den Mittelpunkt (andere qx,qy,qz Berechnungen)

            // Controls group
            "N-;EditBFactor;3|0.01|999;;1",  //??
            "T-;CheckBoxWAXS;;;0",  //??
            "N-;P1;;;0",       //?? wird in formpq verwendet

            // Pixel Manipulation group
            "NF;iso;5|0|1000;Offset for the radial intensity;0",
            "NF;I0;5|0|100000;Multiplier for the radial intensity;10000",
            "NF;Base;5|-10000|10000;Constant background added to the calculated pixel value;0",  // Eigentlich 1e-10
        };
    return slGUI;
}

void SC_Calc_GENERIC::prepareData(_dataGetter dg, bool use1d)
{
    DT( qDebug() << "prepareData()" << use1d );
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
    (*dg)( "HKLmax", val );              calc->setHKLmax( val.value );
    (*dg)( "EditQmax", val );            calc->setQMax( val.value );
    if ( use1d )
    {
        (*dg)( "EditQsteps", val );      calc->setGridPoints( val.value, true );
        (*dg)( "EditQmin", val );        calc->setQMin( val.value );
    }
    else
    {
        (*dg)( "GridPoints", val );      calc->setGridPoints( val.value, false );
        (*dg)( "RadioButtonQ1", val );   calc->setRadQ1( val.value > 0 );   // Weil über InputValues<>
        (*dg)( "RadioButtonQ2", val );   calc->setRadQ2( val.value > 0 );
        (*dg)( "RadioButtonQ4", val );   calc->setRadQ4( val.value > 0 );
        (*dg)( "ExpandImage", val );     calc->setExpandImage( val.value > 0 );
    }

    // Experiment group
    (*dg)( "EditPixelNoX", val );        calc->setpixnox( val.value );
    (*dg)( "EditPixelNoY", val );        calc->setpixnoy( val.value );
    (*dg)( "EditPixelX", val );          calc->setpixx( val.value );    // früher Wert in Meter, jetzt in Millimeter übergeben
    (*dg)( "EditPixelY", val );          calc->setpixy( val.value );    // -"-
    (*dg)( "EditDet", val );             calc->setdet( val.value );
    (*dg)( "EditWavelength", val );      calc->setwave( val.value );
    (*dg)( "BeamPosX", val  );
    (*dg)( "BeamPosY", val1 );           calc->setBeamStop( val.value, val1.value );
    (*dg)( "CenterBeam", val );          calc->setUseBeamStop( val.checked );

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
