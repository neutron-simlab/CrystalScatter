#!/usr/bin/python3.8

# Will not work with python3.6 or less !


# Testroutine to call the CrystalScatter with MCP parameters
# Adapted to the iff1585 linux system (07. Jul. 2025)

# Definition of the CrystalScatter executable
#crystalscatter = "/opt/sas_scatter2/bin/sas_scatter2Cons"
crystalscatter = r"C:\SimLab\sas-crystal\build-sas_scatter2Cons-Desktop_Qt_6_7_2_MinGW_64_bit-Debug\debug\sas_scatter2Cons.exe"

parfile = r"C:\SimLab\CrystalScatter\Mcp\testmcppar.txt"
imgfile = r"C:\SimLab\CrystalScatter\Mcp\testmcpout.png"

import subprocess


# Output from the Chatbot as a list of strings:
keyval = [
    "Experiment_name=Disks",
    "base=0.0",
    "beamposx=57.0",
    "beamposy=31.0",
    "cbinterior=homogeneous",
    "cbparticle=disk",
    "ceff=0.01",
    "ceffcyl=0.0",
    "dbeta=0.4",
    "debyescherrer=False",
    "det=10.0",
    "gridpoints=100.0",
    "i0=1000.0",
    "iso=0.0",
    "length=2.0",
    "ltype=None",
    "ordis=isotropic",
    "peakpar=0.0",
    "phi=0.0",
    "pixelnox=128.0",
    "pixelnoy=128.0",
    "pixelx=1.0",
    "pixely=1.0",
    "qmax=2.0",
    "radius=4.0",
    "rbpara=False",
    "reff=0.0",
    "rotalpha=0.0",
    "rotphi=0.0",
    "rottheta=0.0",
    "sigma=0.1",
    "sigmal=0.1",
    "theta=0.0",
    "ucalpha=90.0",
    "ucbeta=90.0",
    "ucgamma=90.0",
    "ucn1=1.0",
    "ucn2=0.0",
    "ucn3=0.0",
    "ucpsi=0.0",
    "wavelength=0.154"
]

# Convert it in a single string without blanks (use ; as separator)
kvstr = ';'.join(keyval)

# Setup the argument array
args = [ crystalscatter,      # executable
         #"--mcpval", kvstr,  # (inp) list of values from chatbot
         "--mcpinp", parfile, # (inp) filename for the input
         #"--mcpimg64",       # (out) output Base64 image string
         "--mcpimg", imgfile, # (out) filename for the generated output image
         "--threads", "0"     # use the GPU if available (remove to use all cpu cores and no gpu)
         ]
# call the executable with "--help" to get all parameter informations


# Call it
cp = subprocess.run(args, capture_output=True, text=True)
if cp.returncode != 0:
    # This is an error
    print("Returncode:", cp.returncode)
    print(cp.stderr)  # for more explanations
    exit(0)
else:
    # Normal finish...
    print("--- Output ---")
    print(cp.stderr) # Logging to console, I use stderr because there is no buffering
    #      -> 105.381                                    is the runtime in ms
    #     Done in 0 sec =00:00:00
    
    print("--- Image ---")
    tmp = cp.stdout # Base64 Image string output
    print(len(tmp))  # in this example:  26285
    print(tmp)

