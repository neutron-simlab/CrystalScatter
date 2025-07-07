# Testroutine to call the CrystalScatter with MCP parameters

# Current working directory:
#print(os.getcwd())
#C:\SimLab\CrystalScatter\Mcp

# Definition of the CrystalScatter executable
# Final version:
#crystalscatter = "..\\build-sas_scatter2Cons-Desktop_Qt_6_7_2_MinGW_64_bit-Debug\\debug\\sas_scatter2Cons.exe"
# Development:
crystalscatter = "..\\..\\sas-crystal\\build-sas_scatter2Cons-Desktop_Qt_6_7_2_MinGW_64_bit-Debug\\debug\\sas_scatter2Cons.exe"
#crystalscatter =  "C:\\SimLab\\sas-crystal\\build-sas_scatter2Cons-Desktop_Qt_6_7_2_MinGW_64_bit-Debug\\debug\\sas_scatter2Cons.exe"

import subprocess


# Output from the Chatbot:
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

args = [ crystalscatter, "--mcpval", kvstr, "--mcpimg64" ]
#args = [ crystalscatter, "--help" ]

cp = subprocess.run(args, capture_output=True, text=True)
if cp.returncode != 0:
    print("Returncode:", cp.returncode)
    exit(0)
else:
    print("--- Output ---")
    print(cp.stderr) # Logging to console, I use stderr because there is no buffering
    print("--- Image ---")
    print(cp.stdout) # Base64 Image

