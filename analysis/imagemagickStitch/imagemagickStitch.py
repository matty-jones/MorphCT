import sys
import subprocess

if __name__ == "__main__":
    montageDims = sys.argv[1]
    imagesToStitch = sys.argv[2:]
    # First load all the files
    for ID, image in enumerate(imagesToStitch):
        # Load image using supersampling to keep the quality high, and add annotation
        subprocess.call(["convert", "-density", "500", image, "-resize", "20%", "-font", "Arial-Black", "-pointsize", "10", "-gravity", "NorthWest", "-annotate", "0", str(ID+1) + ")", str(ID) + ".png"])
    # Create montage
    montage = subprocess.Popen(["montage", "-mode", "concatenate", "-tile", montageDims, "*.png", "miff:-"], stdout=subprocess.PIPE)
    convert = subprocess.call(["convert", "miff:-", "-density", "500", "-resize", "2000x", "output.png"], stdin=montage.stdout)
    montage.wait()
