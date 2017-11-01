import sys
import os
import subprocess
import glob

def stitchImages(montageDims, imagesToStitch):
    # First delete previous output to prevent ImageMagick issue
    try:
        os.remove("./output.png")
    except:
        pass
    # Then load all the files
    print("Converting images and adding annotations...")
    for ID, image in enumerate(imagesToStitch):
        IDStr = "%04d" % (ID)
        # Load image using supersampling to keep the quality high, and add annotation
        subprocess.call(["convert", "-density", "500", image, "-resize", "20%", "-font", "Arial-Black", "-pointsize", "10", "-gravity", "NorthWest", "-annotate", "0", str(ID+1) + ")", IDStr + "_temp.png"])
    # Create montage
    print("Creating montage...")
    montage = subprocess.Popen(["montage", "-mode", "concatenate", "-tile", montageDims, "*.png", "miff:-"], stdout=subprocess.PIPE)
    print("Exporting montage...")
    convert = subprocess.call(["convert", "miff:-", "-density", "500", "-resize", "2000x", "output.png"], stdin=montage.stdout)
    montage.wait()
    print("Removing temporary files...")
    for fileName in glob.glob("*_temp.png"):
        os.remove(fileName)
    print("Montage created and saved at ./output.png")


if __name__ == "__main__":
    montageDims = sys.argv[1]
    imagesToStitch = sys.argv[2:]
    stitchImages(montageDims, imagesToStitch)
