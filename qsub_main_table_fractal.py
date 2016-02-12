#!/usr/bin/env python
import subprocess


def torque_nageru():
    sketch_k = 128
    pass_num = 1000
    upper_param = 1.0
    graph_names = [
        # "'flower 29526 1 2'",
        # "'flower 88575 1 2'",
        # "'flower 265722 1 2'",

        # "'flower 10924 1 3'",
        # "'flower 43692 1 3'",
        # "'flower 699052 1 3'",

        # "'flower 11720 1 4'",
        # "'flower 58595 1 4'",
        # "'flower 292970 1 4'",

        # "'shm 12501 5 2'",
        # "'shm 62501 5 2'",
        # "'shm 312501 5 2'",

        # "'shm 67229 5 3'",

        # "'flower 43692 2 2'",
        # "'flower 699052 2 2'",
        # "'flower 11720 2 3'",

        # "'flower 58595 2 3'",
        # "'flower 292970 2 3'",

        # "'flower 37326 2 4'",
        # "'flower 223950 2 4'",

        # "'flower 37326 3 3'",
        # "'flower 223950 3 3'",

        # "'flower 14007 3 4'",
        # "'flower 686287 3 4'",


    ]

    methods = [
        "sketch",
        "memb",
        "cbb",
        "burning",
        "coloring",
    ]

    exp_tag = "main_table"
    rad_max = 1000000

    # # MEMB
    for method in methods:
        for graph_name in graph_names:
            command = "/home/kenkoooo/fractal-dimension/agl/bin/box_cover --force_undirected"
            command = command + " --type gen "
            command = command + " --graph " + graph_name
            command = command + " --rad_analytical "
            command = command + " --rad_max " + str(rad_max)
            command = command + " --method " + method
            command = command + " --exp_tag " + exp_tag
            command = command + " --sketch_k " + str(sketch_k)
            command = command + " --pass " + str(pass_num)
            command = command + " --upper_param " + str(upper_param)
            command = command + " --final_coverage " + str(1.0)

            job_name = exp_tag + "-" + method + "-" + \
                graph_name.replace("'", "").replace(" ", "-")
            p1 = subprocess.Popen(["echo", command], stdout=subprocess.PIPE)
            p2 = subprocess.Popen(
                ["qsub", "-l", "walltime=24:00:00,nodes=node06:ppn=2", "-N", job_name], stdin=p1.stdout)
            p1.stdout.close()
            output = p2.communicate()[0]

if __name__ == "__main__":
    torque_nageru()
