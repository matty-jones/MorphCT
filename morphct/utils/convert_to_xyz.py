import sys


def convert_inp(file_name):
    with open(file_name, "r") as inp_file:
        inp_data = inp_file.readlines()
    xyz_data = []
    for line in inp_data:
        if line[0] == " ":
            xyz_data.append(line)
    xyz_data.insert(0, "{:d}\n".format(len(xyz_data)))
    xyz_data.insert(1, "comment_line\n")
    with open(file_name.replace(".inp", ".xyz"), "w+") as xyz_file:
        xyz_file.writelines(xyz_data)


def convert_out(file_name):
    with open(file_name, "r") as out_file:
        out_data = out_file.readlines()
    xyz_data = []
    for line in out_data:
        if line[0] == "|":
            xyz_data.append(line)
    # Remove the header lines
    while "*" not in xyz_data[0]:
        xyz_data.pop(0)
    xyz_data.pop(0)
    # Remove final lines and ORCA formatting
    pop_list = []
    for line_no, line in enumerate(xyz_data):
        if ("*" in line) or ("end" in line) or (len(line.split("> ")[-1][:-1]) < 3):
            pop_list.append(line_no)
        else:
            xyz_data[line_no] = line.split("> ")[-1]
    for line_no in sorted(pop_list, reverse=True):
        xyz_data.pop(line_no)
    xyz_data.insert(0, "{:d}\n".format(len(xyz_data)))
    xyz_data.insert(1, "comment_line\n")
    with open(file_name.replace(".out", ".xyz"), "w+") as xyz_file:
        xyz_file.writelines(xyz_data)


def main():
    list_of_files = sys.argv[1:]
    if len(list_of_files) < 1:
        print("No files requested to convert!")
        exit()
    for file_name in list_of_files:
        print("Converting {:s}...".format(file_name))
        if ".inp" in file_name:
            convert_inp(file_name)
        elif ".out" in file_name:
            convert_out(file_name)
        else:
            print("Unknown file type for {:s}, skipping...".format(file_name))
    print("Conversion tasks completed!")


if __name__ == "__main__":
    main()
