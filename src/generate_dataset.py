import os
import sys
import argparse
import shutil
import filecmp
import subprocess
import logging

logging.basicConfig(
    format="[%(filename)s:%(lineno)s - %(funcName)20s() ] %(message)s",
    stream=sys.stderr, level=logging.DEBUG)

# Assumptions
# 1. this file will be executed from within the src folder
# 2. pwd with the esso/src


def copy_file_if_not_present_or_older(base_path, new_path):

    try:
        # if the file is not present then copy it
        if not os.path.exists(new_path):
            shutil.copy(base_path, new_path)
        # file exists but it is not same as the base path, then replace
        elif not filecmp.cmp(base_path, new_path):
            shutil.copy(base_path, new_path)
        else:
            pass
    except:
        logging.error('Failed to copy file for {} {}'.format(
            base_path, new_path))


def execute(cmd):
    """
    execute a command and provide an iterable to read its output
    :param cmd:
    :return:
    """
    logging.debug('params {}'.format(" ".join(cmd)))
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                             universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line
    popen.stdout.close()
    return_code = popen.wait()
    print(return_code)
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)


def generate_dataset(as_num, arrival_rate, sfc_lifetime, replace=False):
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
    logging.debug('params {} {} {} {}'.format(as_num, arrival_rate,
                                              sfc_lifetime, replace))
    dataset_dir_name = as_num + '_' + arrival_rate + \
                       '_' + sfc_lifetime
    dataset_dir_path = '../data/' + dataset_dir_name
    try:
        if os.path.isdir(dataset_dir_path) and not replace:
            return
    except:
        logging.error('Failed to check dataset dir existence for {} {} {} {}'.
            format(as_num, arrival_rate, sfc_lifetime, replace))

    #############################
    # now moving onto the tasks #
    #############################

    # create the directory if not present
    try:
        if not os.path.isdir(dataset_dir_path):
            os.mkdir(dataset_dir_path)
    except:
        logging.error('Failed to check or create dataset dir for {} {} {} {}'.
                      format(as_num, arrival_rate, sfc_lifetime, replace))

    # check for vnf_type file
    vnf_type_file_path = os.path.join(dataset_dir_path, 'vnf_types.dat')
    copy_file_if_not_present_or_older('../data/set0/vnf_types.dat',
                                      vnf_type_file_path)

    # check for co_topology.dat file
    co_topology_file_path = os.path.join(dataset_dir_path, 'co_topology.dat')
    copy_file_if_not_present_or_older('../data/rocketfuel_v2/'+
                                      as_num+'.co', co_topology_file_path)

    # run process_topology.o
    try:
        for output in execute(["./process_topology.o", dataset_dir_path]):
            pass
    except:
        logging.error('Failed to process topology for {} {} {} {}'.format(
            as_num, arrival_rate, sfc_lifetime, replace))


    # run traffic_generator
    try:
        current_dir = os.getcwd()
        os.chdir('../data/rocketfuel/')
    except:
        logging.error('Failed to store and change dir for {} {} {} {}'.
                      format(as_num, arrival_rate, sfc_lifetime, replace))

    try:
        for output in execute(['python', 'traffic_generator.py',
                           '--as_no', as_num,
                           '--arrival_rate', arrival_rate,
                           '--sfc_lifetime', sfc_lifetime]):
            pass
    except:
        logging.error('Failed to generate traffic for {} {} {} {}'.format(
            as_num, arrival_rate, sfc_lifetime, replace))

    try:
        os.chdir(current_dir)
    except:
        logging.error('Failed to chage dir back to src for {} {} {} {}'.
                      format(as_num, arrival_rate, sfc_lifetime, replace))

    # copy timeslot file
    timeslots_file_path = os.path.join(dataset_dir_path, 'timeslots.dat')
    copy_file_if_not_present_or_older('../data/rocketfuel/'+
                                      as_num+'.timeslots',
                                      timeslots_file_path)


def get_args():
    parser = argparse.ArgumentParser(formatter_class=
                                     argparse.ArgumentDefaultsHelpFormatter)
    # if the `-r` option is provide then any existing directory under
    # the `runs` folder with the same run id will be replaced
    parser.add_argument('-r', '--replace', action='store_true',
                        help="replace existing dataset dir")
    parser.add_argument("--as_num", type = str, required = True)

    return parser.parse_args()


if __name__ == "__main__":

    #################################

    # arrival rates will be divided by 100
    sfc_arrival_rate_start = 10
    sfc_arrival_rate_increment = 1
    sfc_arrival_rate_datapoints = 10 # 10? how many data points
    sfc_arrival_rate_end = sfc_arrival_rate_start + \
                           sfc_arrival_rate_increment * \
                           sfc_arrival_rate_datapoints

    sfc_lifetime_start = 3600
    sfc_lifetime_increment = 1800
    sfc_lifetime_datapoints = 20 # 20?
    sfc_lifetime_end = sfc_lifetime_start + \
                       sfc_lifetime_increment * sfc_lifetime_datapoints

    #################################

    args = get_args()
    for sfc_arrival_rate in range(sfc_arrival_rate_start,
                                  sfc_arrival_rate_end,
                                  sfc_arrival_rate_increment):
        for sfc_lifetime in range(sfc_lifetime_start,
                                  sfc_lifetime_end,
                                  sfc_lifetime_increment):
            try:
                generate_dataset(args.as_num,
                                 '{:4.2f}'.format(sfc_arrival_rate * 0.01),
                                 str(sfc_lifetime), args.replace)
            except Exception as e:
                logging.error('Failed for {} {}'.format(
                    sfc_arrival_rate*0.01, sfc_lifetime))
                logging.error(str(e))

