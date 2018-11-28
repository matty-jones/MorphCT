import os
import subprocess
import argparse


def get_parameter_lines(morphology_name, input_xml):
    with open("par.py", "r") as parameter_file:
        lines = parameter_file.readlines()
        for line_number, line in enumerate(lines):
            if "<INPUTXML>" in line:
                lines[line_number] = "".join(["'", input_xml, "'"]).join(
                    line.split("<INPUTXML>")
                )
    new_parameter_location = os.path.join(
        ".", "submissions", "".join(["par_", morphology_name, ".py"])
    )
    with open(new_parameter_location, "w+") as new_parameter_file:
        new_parameter_file.writelines(lines)
    return new_parameter_location


def get_slurm_sub_lines(morphology_name, input_xml, param_file_loc):
    with open("submit.sh", "r") as slurm_file:
        lines = slurm_file.readlines()
        for line_number, line in enumerate(lines):
            if "<INPUTXML>" in line:
                lines[line_number] = morphology_name.join(line.split("<INPUTXML>"))
            elif "<INPUTPAR>" in line:
                lines[line_number] = param_file_loc.join(line.split("<INPUTPAR>"))
    new_submission_location = os.path.join(
        ".", "submissions", "".join(["sub_", morphology_name, ".sh"])
    )
    with open(new_submission_location, "w+") as new_submission_file:
        new_submission_file.writelines(lines)
    return new_submission_location


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    args, files_list = parser.parse_known_args()
    if len(files_list) == 0:
        files_list = os.listdir("./inputs")
    else:
        files_list = [os.path.split(file_name)[-1] for file_name in files_list]
    for input_xml in files_list:
        morphology_name = os.path.splitext(input_xml)[0]
        print("Writing parameter file and submission script for", morphology_name)
        parameter_location = get_parameter_lines(morphology_name, input_xml)
        new_submission_location = get_slurm_sub_lines(
            morphology_name, input_xml, parameter_location
        )
        subprocess.call(["sbatch", new_submission_location])
