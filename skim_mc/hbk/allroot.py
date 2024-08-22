import os
for root, dirs, files in os.walk("."):  
    for filename in files:
        if ".hist" in filename:
            os.system("h2root {}".format(filename))
os.system(r"hadd sup.root *.root")
os.system(r"mv sup.root ~/Mat_belle/root")
os.system(r"rm -f *.root")
