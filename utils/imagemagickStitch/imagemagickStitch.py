import sys
import os
import subprocess
import glob
import argparse

def stitchImages(montageDims, imagesToStitch, morphologyName, title):
    # First delete previous output to prevent ImageMagick issue
    if title is None:
        title = morphologyName
    try:
        os.remove("./" + title.replace(" ", "_") + ".png")
    except:
        pass
    # Then load all the files
    print("Converting images and adding annotations...")
    for ID, image in enumerate(imagesToStitch):
        IDStr = "%04d" % (ID)
        if image[-4:] == '.pdf':
            ## Load image using supersampling to keep the quality high, and "crop" to add a section of whitespace to the left
            #subprocess.call(["convert", "-density", "500", image, "-resize", "20%", "-gravity", "West", "-bordercolor", "white", "-border", "7%x0", IDStr + "_crop.png"])
            ## Now add the annotation
            #subprocess.call(["convert", IDStr + "_crop.png", "-font", "Arial-Black", "-pointsize", "72", "-gravity", "NorthWest", "-annotate", "0", str(ID+1) + ")", IDStr + "_temp.png"])
            subprocess.call(["convert", "-density", "500", image, "-resize", "20%", "-font", "Arial-Black", "-pointsize", "72", "-gravity", "NorthWest", "-splice", "0x10%", "-page", "+0+0", "-annotate", "0", str(ID+1) + ")", IDStr + "_temp.png"])
        else:
            # Rasterized format, so no supersampling
            #subprocess.call(["convert", image, "-gravity", "West", "-bordercolor", "white", "-border", "7%x0", IDStr + "_crop.png"])
            # Now add the annotation
            subprocess.call(['convert', image, '-font', 'Arial-Black', '-pointsize', '72', '-gravity', 'NorthWest', '-splice', '0x10%', '-page', '+0+0', '-annotate', '0', str(ID+1) + ')', IDStr + '_temp.png'])
    # Create montage
    print("Creating montage...")
    montage = subprocess.Popen(["montage", "-mode", "concatenate", "-tile", montageDims, "*_temp.png", "miff:-"], stdout=subprocess.PIPE)
    print("Exporting montage...")
    convert = subprocess.call(["convert", "miff:-", "-density", "500", "-resize", "2000x", "-font", "Arial-Black", "-pointsize", "10", "-gravity", "North", "-bordercolor", "white", "-border", "0x100", "-annotate", "0", title, title.replace(" ", "_") + ".png"], stdin=montage.stdout)
    montage.wait()
    print("Removing temporary files...")
    for fileName in glob.glob("*_temp.png") + glob.glob("*_crop.png"):
        os.remove(fileName)
    print("Montage created and saved at ./" + title.replace(" ", "_") + ".png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--dimensions", default="x1", required=False, help="The dimensions flag is used to specify the montage dimensions. The format should be '2x3', which corresponds to a montage with 2 columns and 3 rows. Dimensions can be omitted such as '3x', which will create a montage with 3 columns and as many rows as required based on the number of input figures. Default is a single row of all input images.")
    parser.add_argument("-t", "--title", default=None, required=False, help="Set a custom title for the montage. If not set, will be assigned based on the enclosing directory.")
    args, imagesToStitch = parser.parse_known_args()
    currentDir = os.getcwd()
    morphologyName = currentDir[currentDir.rfind("/") + 1:]
    stitchImages(args.dimensions, imagesToStitch, morphologyName, args.title)
