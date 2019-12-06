# Setup.py A small script to copy the files needed for anaylsis testing, but not for development, they are stored in the dir above
# Ben Lorentz 11.20.2019
import os

whereAmI = os.getcwd()
os.chdir("..")
os.system("cp -a sample/. LorentzThesis19/")
os.chdir(whereAmI)
print("All set")
