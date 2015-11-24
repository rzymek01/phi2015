#!/usr/bin/python
import subprocess
import os

print os.getcwd()

OUTPUT_PATH = "/macierz/home/137396rm/cuda/msg-pass-sim-tests/"
INPUT_PATH = "/macierz/home/137396rm/cuda/msg-pass-sim-tests/"
PROGRAM_PATH = "/macierz/home/137396rm/cuda/msg-pass-sim/Release/msg"


def execute(specs_id, v, t, d, test):
    cmd = PROGRAM_PATH + str(v) + " " + str(t) + " " + str(d) + " <" + INPUT_PATH + "test" + str(test) + ".in >"\
        + OUTPUT_PATH + "test_" + specs_id + "_" + str(v) + "_" + str(d) + "_" + str(t) + "_" + str(test) + ".out 2>"\
        + OUTPUT_PATH + "test_" + specs_id + "_" + str(v) + "_" + str(d) + "_" + str(t) + "_" + str(test) + ".time"

    subprocess.call(cmd, shell=True)
#    print(cmd)


def run(specs):
    for i in specs:
        # print(i)
        specs_id = specs[i]["id"]
        specs_name = specs[i]["name"]
        specs_cmd = specs[i]["cmd"]

        for t in specs_cmd["t"]:
            # print(t)
            for d in specs_cmd["d"]:
                # print(d)
                for v in specs_cmd["v"]:
                    # print(v)
                    for test in specs_cmd["test"]:
                        # print(test)
                        execute(specs_id, v, t, d, test)


def find_diff(specs, v1, v2):
    for i in specs:
        # print(i)
        specs_id = specs[i]["id"]
        specs_name = specs[i]["name"]
        specs_cmd = specs[i]["cmd"]
        specs_v = specs[i]["cmd"]["v"]

        specs_v_len = len(specs_v)
        if specs_v_len == 1:
            # print "nothing to do in diff"
            return
        else:
            # print "check diff"
            for t in specs_cmd["t"]:
                # print(t)
                for d in specs_cmd["d"]:
                    # print(d)
                    for test in specs_cmd["test"]:
                        # print(test)
                        path = OUTPUT_PATH
                        diff_cmd = "diff -q " + path + "test_"\
                                   + specs_id + "_" + str(v1) + "_" + str(d) + "_" + str(t) + "_" + str(test)\
                                   + ".out "\
                                   + path + "test_"\
                                   + specs_id + "_" + str(v2) + "_" + str(d) + "_" + str(t) + "_" + str(test)\
                                   + ".out"
                        try:
                            diff_output = subprocess.check_output(diff_cmd, shell=True)
                        except subprocess.CalledProcessError as e:
                            print e.output
                            exit("Outputs are different. Something wrong with program.")

# main
if __name__ == "__main__":
    most_optimized_version = 4
    faster_device = 0
    threads_per_block = 128

    specs = {
        # "A": {
        #     "id": "A3",
        #     "name": "Test A. Execution time vs number of threads per block",
        #     "cmd": {
        #         "t": (64, 128, 256, 512, 1024),
        #         "d": (faster_device,),
        #         "v": (most_optimized_version,),
        #         "test": (1000, '1000000b')
        #     }
        # },
        # "B": {
        #     "id": "B3",
        #     "name": "Test B. Execution time vs optimizations",
        #     "cmd": {
        #         "t": (threads_per_block,),
        #         "d": (faster_device,),
        #         "v": (1, 2, 3, most_optimized_version),
        #         "test": (100, 1000, '100b', '1000b', '10000b', '100000b', '1000000b')
        #     }
        # },
        # "C": {
        #     "id": "C3",
        #     "name": "Test C. Execution time vs size of input data (number of vertices and edges)",
        #     "cmd": {
        #         "t": (threads_per_block,),
        #         "d": (faster_device,),
        #         "v": (most_optimized_version,),
        #         "test": (10, 100, 1000, '10b', '100b', '1000b', '10000b', '100000b', '1000000b')
        #     }
        # },
        # "D": {
        #     "id": "D3",
        #     "name": "Test D. Execution time vs device",
        #     "cmd": {
        #         "t": (threads_per_block,),
        #         "d": (0, 1),
        #         "v": (most_optimized_version,),
        #         "test": (100, 1000, '100000b', '1000000b')
        #     }
        # },
        # "E": {
        #     "id": "E3",
        #     "name": "Test E. Execution time vs size of input data, number of threads per block and optimizations",
        #     "cmd": {
        #         "t": (64, 128, 256, 512, 1024),
        #         "d": (faster_device,),
        #         "v": (1, most_optimized_version),
        #         "test": ('100', '1000', '100b', '1000b', '10000b', '100000b', '1000000b')
        #     }
        # },
        "F": {
            "id": "F3",
            "name": "Test F. Execution time vs size of input data, number of threads per block and device",
            "cmd": {
                "t": (64, 128, 256, 512, 1024),
                "d": (0, 1),
                "v": (most_optimized_version,),
                "test": ('100', '1000', '100b', '1000b', '10000b', '100000b', '1000000b')
            }
        },
    }

    run(specs)
    # find_diff(specs, 1, most_optimized_version)
