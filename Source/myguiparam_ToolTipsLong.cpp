/**
 * Program: sas_scatter2
 *  Author: Michael Wagener, FZJ, JCNS-1
 *          Stephan Förster, FZJ, JCNS-1 wrote the original Scatter in Pascal.
 *    Date: March 2023
 *    File: sc_mainToolTipsLong.cpp
 *          Contains a longer version for the tool tips of all labels.
 *          These texts are taken from the Chatbot system
 *          (https://vipr-staging.chatbot.hw.ag/demo/a8669488-7e02-4bc8-ab03-cc39fba84d41/fullscreen)
 *          (go.fzj.de/crystal-scatter-chatbot)
 */


#include "myguiparam.h"
#include <QDebug>
#include <QRegularExpression>


bool myGuiParam::localSetLongerTooltip(QWidget *w)
{
    static QString ltt[] = // LongerToolTip
        {
            // Lattice group
            "LType", "The 'LType' parameter specifies the lattice type in the Crystal Scatter simulation, which is crucial"
                     " for defining the atomic arrangement within a crystal. It influences the scattering patterns calculated"
                     " by the software, as different lattice types (e.g., BCC, FCC, HCP) have distinct symmetries and atomic"
                     " positions. Selecting the appropriate lattice type is essential for accurate simulations, ensuring that"
                     " the results align with the actual material being analyzed. Typical values for 'ltype' include various"
                     " arrangements such as Lamellae, Gyroid, and Orthorhombic, reflecting the diverse structures found in"
                     " crystalline materials.",
            "LType1d", "",

            "uca", "",
            "ucb", "",
            "ucc", "",
            "ucalpha", "",
            "ucbeta", "",
            "ucgamma", "",
            "EditCeff", "",
            "CheckBoxTwinned", "",
            "EditCeffcyl", "",
            "reff", "",
            "acpl", "",
            "bcpl", "",
            "ifluc", "",
            "rfluc", "",

            // Particle group
            "ComboBoxParticle", "",
            "Ordis", "",
            "ComboBoxInterior", "",

            //(* Particle Interior *) - 1D
            "ComboBoxInterior1d", "",

            "EditRadius", "",
            "EditRadiusi", "",
            "EditSigma", "",
            "EditDbeta", "The parameter 'dbeta' denotes the average scattering angle, currently set at 20.0° with a range"
                         " from 0° to 360°. In crystallography, this angle is pivotal as it defines the direction in which"
                         " scattered radiation is measured, directly influencing the resulting diffraction pattern. A smaller"
                         " angle, like 20.0°, typically corresponds to the observation of lower-order reflections, which can"
                         " provide insights into the overall symmetry and periodicity of the crystal structure. The choice of"
                         " 'dbeta' is crucial for optimizing experimental setups, as it affects the resolution and detail of"
                         " the data collected. Scattering patterns at various angles can reveal different aspects of the"
                         " crystal's structure, making this parameter essential for accurate structural characterization.",
            "Length", "",
            "SigmaL", "",
            "Alpha", "",
            "EditRho", "",
            "RotAlpha", "",

            // Peak Shape group
            "ComboBoxPeak", "",
            "EditPeakPar", "",
            "EditDebyeWaller", "",
            "EditAzi", "",
            "EditDomainSize", "",
            "RadButDebyeScherrer", "",
            "RadioButtonPara", "",
            "VAx1", "",
            "VAx2", "",
            "VAx3", "",
            "Ay1", "",
            "Ay2", "",
            "Ay3", "",
            "Az1", "",
            "Az2", "",
            "Az3", "",
            "SigX", "",
            "SigY", "",
            "SigZ", "",

            // Orientation group
            "ucpsi", "",
            "ucn1", "",
            "ucn2", "",
            "ucn3", "",
            "theta", "",
            "phi", "",
            "rotTheta", "",
            "rotPhi", "",

            // Calculation group
            "GridPoints", "",
            "HKLmax", "",
            "EditQmax", "",
            "CalcQmax", "",
            "EditQmaxData", "",
            "EditQmaxPreset", "",
            "EditQmin", "",
            "EditQsteps", "",
            "Editqz", "",

            // Experiment group
            "EditPixelNoX", "",
            "EditPixelNoY", "",
            "EditPixelX", "",
            "EditPixelY", "",
            "EditDet", "The 'EditDet' parameter represents the distance from the sample to the detector in a scattering experiment."
                       " It is crucial for determining the resolution of the scattering data. A distance of 10.0 meters is typical"
                       " and serves as the default value, with a permissible range from 0.001 to 100 meters. Adjusting this distance"
                       " affects the angular resolution and intensity of the scattering pattern, making it vital for optimizing"
                       " experimental setups in crystallography.",
            "EditWavelength", "",
            "BeamPosX", "",
            "BeamPosY", "",
            "CenterBeam", "",
            "CenterMidpoint", "",

            // Controls group
            "EditBFactor", "",
            "CheckBoxWAXS", "",
            "P1", "",

            // Pixel Manipulation group
            "iso", "",
            "I0", "The parameter 'I0' is a multiplier for the radial intensity in scattering simulations. It is set at a current"
                  " and default value of 1000.0, with a range from 0 to 100,000. In crystallography, 'i0' affects the overall"
                  " intensity of the scattering pattern, allowing researchers to enhance the visibility of features in the data."
                  " Adjusting 'i0' can optimize experimental setups, improving signal detection and structural characterization.",
            "Base", "",

            "", "" // EndOfList
        };

    int i=0;
    QString cmp = w->objectName().mid(3);
    while ( ltt[i].length() > 0 )
    {
        if ( ltt[i] == cmp )
        {
            if ( ltt[i+1].isEmpty() ) return false; // Noch nicht definiert
            w->setToolTip(ltt[i+1]);
            return true;
        }
        i+=2;
    }
    qDebug() << "localSetLongerTooltip - not found" << cmp;
    return false;
}
