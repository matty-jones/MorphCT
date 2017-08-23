from hoomd_script import *
import sys

if __name__ == "__main__":
    fileName = sys.argv[1]
    system = init.read_xml(filename = fileName)
    dump.dcd(filename = fileName.replace(".xml", ".dcd"), overwrite = True, period = 1)
    run(1)
    print("DCD file created for", fileName, "as", fileName.replace(".xml", ".dcd"))
