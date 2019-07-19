import os
import argparse
import shutil
import filecmp

# Assumptions
# 1. this file will be executed from within the src folder
# 2. pwd with the esso/src


def copy_file_if_not_present_or_older(base_path, new_path):

    # if the file is not present then copy it
    if not os.path.exists(new_path):
        shutil.copy(base_path, new_path)
    # file exists but it is not same as the base path, then replace
    elif not filecmp.cmp(base_path, new_path):
        shutil.copy(base_path, new_path)
    else:
        pass


def prepare_dataset_dir(as_num, arrival_rate, sfc_lifetime, replace=False):
    """
    check whether dataset dir exists.
    If exists and `-r` option not provided then exit.

    If the dir does not exist or `-r` option provide then create the dir.
    Perform the following tasks:
    1. check and update co topology file if needed
    2. check and update vnf_type file if needed
    3. run process_topology
    4. generate timeslots.dat file
    :return:
    """
    dataset_dir_name = as_num + '_' + arrival_rate + '_' + sfc_lifetime
    dataset_dir_path = '../data/' + dataset_dir_name
    if os.path.isdir(dataset_dir_path) and not replace:
        return

    #############################
    # now moving onto the tasks #
    #############################

    # create the directory if not present
    if not os.path.isdir(dataset_dir_path):
        os.mkdir(dataset_dir_path)

    # check for vnf_type file
    vnf_type_file_path = os.path.join(dataset_dir_path, 'vnf_types.dat')
    copy_file_if_not_present_or_older('../data/set0/vnf_types.dat', vnf_type_file_path)

    # check for co_topology.dat file
    co_topology_file_path = os.path.join(dataset_dir_path, 'co_topology.dat')
    copy_file_if_not_present_or_older(as_num+'.co', co_topology_file_path)

    # run process_topology.o


def generate_dataset():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # if the `-r` option is provide then any existing directory under
    # the `runs` folder with the same run id will be replaced
    parser.add_argument('-r', '--replace', action='store_true',
                        help="replace existing dataset dir")
    parser.add_argument("--as_num", type = str, required = True)
    parser.add_argument('--arrival_rate', type=str, default=0.003) # 10 sfc per hour
    parser.add_argument('--sfc_lifetime', type=str, default=5400) # 1 hour 30 minute

if __name__ == "__main__":
    prepare_dataset_dir()

    generate_dataset()