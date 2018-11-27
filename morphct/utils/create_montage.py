import sys
import os
import subprocess as sp
import glob
import argparse


def stitch_images(montage_dims, images_to_stitch, morphology_name, title, save_file):
    # First delete previous output to prevent ImageMagick issue
    if title is None:
        title = morphology_name
    try:
        os.remove("".join("./", title.replace(" ", "_") + ".png"))
    except:
        pass
    # Then load all the files
    print("Converting images and adding annotations...")
    directory = "/".join(images_to_stitch[0].split("/")[:-2])
    for ID, image in enumerate(images_to_stitch):
        ID_str = "{:04d}".format(ID)
        if image[-4:] == ".pdf":
            print("Vectorized input image detected, using supersampling...")
            sp.call(
                [
                    "convert",
                    "-density",
                    "500",
                    image,
                    "-resize",
                    "20%",
                    "-font",
                    "Arial-Black",
                    "-pointsize",
                    "10",
                    "-gravity",
                    "NorthWest",
                    "-splice",
                    "0x10%",
                    "-page",
                    "+0+0",
                    "-annotate",
                    "0",
                    "{:d})".format(ID + 1),
                    os.path.join(directory, "".join([ID_str, "_temp.png"])),
                ]
            )
        else:
            sp.call(
                [
                    "convert",
                    image,
                    "-resize",
                    "500x500",
                    "-font",
                    "Arial-Black",
                    "-pointsize",
                    "72",
                    "-gravity",
                    "NorthWest",
                    "-splice",
                    "100x120",
                    "-page",
                    "+0+0",
                    "-annotate",
                    "+40+0",
                    "{:d})".format(ID + 1),
                    os.path.join(directory, "".join([ID_str, "_temp.png"])),
                ]
            )
    # Create montage
    print("Creating montage...")
    montage = sp.Popen(
        [
            "montage",
            "-mode",
            "concatenate",
            "-tile",
            montage_dims,
            os.path.join(directory, "*_temp.png"),
            "miff:-",
        ],
        stdout=sp.PIPE,
    )
    print("Exporting montage...")
    if save_file is None:
        save_file = os.path.join(directory, "".join([title.replace(" ", "_"), ".png"]))
    convert = sp.call(
        [
            "convert",
            "miff:-",
            "-density",
            "500",
            "-resize",
            "2000x",
            "-font",
            "Arial-Black",
            "-pointsize",
            "10",
            "-gravity",
            "North",
            "-bordercolor",
            "white",
            "-border",
            "0x100",
            "-annotate",
            "0",
            title,
            save_file,
        ],
        stdin=montage.stdout,
    )
    montage.wait()
    print("Removing temporary files...")
    for file_name in glob.glob(os.path.join(directory, "*_temp.png")) + glob.glob(
        os.path.join(directory, "*_crop.png")
    ):
        os.remove(file_name)
    print("Montage created and saved at ", save_file)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d",
        "--dimensions",
        default="x1",
        required=False,
        help=(
            "The dimensions flag is used to specify the montage"
            " dimensions. The format should be '2x3', which corresponds"
            " to a montage with 2 columns and 3 rows. Dimensions can be"
            " omitted such as '3x', which will create a montage with 3"
            " columns and as many rows as required based on the number"
            " of input figures. Default is a single row of all input images."
        ),
    )
    parser.add_argument(
        "-t",
        "--title",
        default=None,
        required=False,
        help=(
            "Set a custom title for the montage. If not set, will be"
            "assigned based on the enclosing directory."
        ),
    )
    parser.add_argument(
        "-s",
        "--save_as",
        default=None,
        required=False,
        help=(
            "Location and name of the output montage. Default is based on"
            "the title in the first directory above the morphology directories"
            "(might be ../ from cwd)."
        ),
    )
    parser.add_argument(
        "-f",
        "--files",
        action="store_true",
        required=False,
        help=(
            "Operate create_montage in file mode, which does not look for"
            "the standard MorphCT operating file structure."
        ),
    )
    args, directories = parser.parse_known_args()
    if args.files:
        images_to_stitch = [
            os.path.join(os.getcwd(), directory) for directory in directories
        ]
        stitch_images(
            args.dimensions, images_to_stitch, "Montage", args.title, args.save_as
        )
    else:
        for directory in directories:
            morphology_name = os.path.split(directory)[1]
            try:
                images_to_stitch = [
                    os.path.join(directory, "figures", figure)
                    for figure in os.listdir(os.path.join(directory, "figures"))
                ]
                if len(images_to_stitch) == 0:
                    raise FileNotFoundError
            except FileNotFoundError:
                continue
            stitch_images(
                args.dimensions,
                images_to_stitch,
                morphology_name,
                args.title,
                args.save_as,
            )


if __name__ == "__main__":
    main()
