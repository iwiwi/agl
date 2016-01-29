#!/usr/bin/env python
import subprocess


def torque_nageru():
    sketch_k = 128
    pass_num = 1000
    upper_param = 1.0
    graph_names = [
        # "'flower 4 2 2'",
        # "'flower 12 2 2'",
        # "'flower 44 2 2'",
        # "'flower 172 2 2'",
        # "'flower 684 2 2'",
        # "'flower 2732 2 2'",
        # "'flower 10924 2 2'",
        # "'flower 43692 2 2'",
        # "'flower 174764 2 2'",
        "'flower 699052 2 2'",
        # "'flower 2796204 2 2'",
        # "'flower 11184812 2 2'",
    ]

    methods = [
        # "sketch",
        # "coloring",
        "memb",
        # "burning",
    ]

    exp_tag = "96G-scalability"
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

            job_name = exp_tag + "-" + method + "-" + \
                graph_name.replace("'", "").replace(" ", "-")
            p1 = subprocess.Popen(["echo", command], stdout=subprocess.PIPE)
            p2 = subprocess.Popen(
                ["qsub", "-l", "walltime=24:00:00,nodes=1:ppn=24", "-N", job_name], stdin=p1.stdout)
            p1.stdout.close()
            output = p2.communicate()[0]

    # Analytical
    # for graph_name in graph_names:
    #     command = "/home/kenkoooo/fractal-dimension/agl/bin/box_cover --force_undirected --type gen"
    #     command = command + " --graph " + graph_name
    #     command = command + " --method " + "analytical"
    #     # method_name = "memb"
    #     # command = command + " --method " + method_name
    #     # job_name = method_name + "-" + graph_name.replace("'", "").replace(" ", "-")
    #     job_name = "analytical-" + \
    #         graph_name.replace("'", "").replace(" ", "-")
    #     p1 = subprocess.Popen(["echo", command], stdout=subprocess.PIPE)
    #     p2 = subprocess.Popen(
    #         ["qsub", "-l", "walltime=240:00:00", "-N", job_name], stdin=p1.stdout)
    #     p1.stdout.close()
    #     output = p2.communicate()[0]

if __name__ == "__main__":
    torque_nageru()
