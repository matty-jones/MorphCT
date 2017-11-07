import sys
import os
import subprocess
import glob

def stitchImages(montageDims, imagesToStitch, morphologyName):
    # First delete previous output to prevent ImageMagick issue
    try:
        os.remove("./output.png")
    except:
        pass
    # Then load all the files
    print("Converting images and adding annotations...")
    for ID, image in enumerate(imagesToStitch):
        IDStr = "%04d" % (ID)
        # Load image using supersampling to keep the quality high, and "crop" to add a section of whitespace to the left
        subprocess.call(["convert", "-density", "500", image, "-resize", "20%", "-gravity", "West", "-bordercolor", "white", "-border", "7%x0", IDStr + "_crop.png"])
        # Load image using supersampling to keep the quality high, and add annotation
        #subprocess.call(["convert", "-density", "500", IDStr + "_crop.png", "-resize", "20%", "-font", "Arial-Black", "-pointsize", "10", "-gravity", "NorthWest", "-annotate", "0", str(ID+1) + ")", IDStr + "_temp.png"])
        subprocess.call(["convert", IDStr + "_crop.png", "-font", "Arial-Black", "-pointsize", "72", "-gravity", "NorthWest", "-annotate", "0", str(ID+1) + ")", IDStr + "_temp.png"])
    # Create montage
    print("Creating montage...")
    montage = subprocess.Popen(["montage", "-mode", "concatenate", "-tile", montageDims, "*_temp.png", "miff:-"], stdout=subprocess.PIPE)
    print("Exporting montage...")
    convert = subprocess.call(["convert", "miff:-", "-density", "500", "-resize", "2000x", "-font", "Arial-Black", "-pointsize", "10", "-gravity", "North", "-bordercolor", "white", "-border", "0x100", "-annotate", "0", morphologyName, "output.png"], stdin=montage.stdout)
    montage.wait()
    print("Removing temporary files...")
    for fileName in glob.glob("*_temp.png") + glob.glob("*_crop.png"):
        os.remove(fileName)
    print("Montage created and saved at ./output.png")


if __name__ == "__main__":
    montageDims = sys.argv[1]
    imagesToStitch = sys.argv[2:]
    currentDir = os.getcwd()
    morphologyName = currentDir[currentDir.rfind("/") + 1:]
    stitchImages(montageDims, imagesToStitch, morphologyName)
