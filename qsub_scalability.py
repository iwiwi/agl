#!/usr/bin/env python
import subprocess


def torque_nageru():
    sketch_k = 128
    pass_num = 1000
    upper_param = 1.0
    graph_names = [
        # "'ba 128 2'",
        # "'ba 256 2'",
        # "'ba 512 2'",
        # "'ba 1024 2'",
        # "'ba 2048 2'",
        # "'ba 4096 2'",
        # "'ba 8192 2'",
        # "'ba 16384 2'",
        # "'ba 32768 2'",
        # "'ba 65536 2'",
        # "'ba 131072 2'",
        # "'ba 262144 2'",
        # "'ba 524288 2'",
        # "'ba 1048576 2'",
        # "'ba 2097152 2'",
        # "'ba 4194304 2'",
        # "'ba 8388608 2'",

        # "'ba 125 2'",
        # "'ba 250 2'",
        # "'ba 260 2'",
        # "'ba 270 2'",
        # "'ba 280 2'",
        # "'ba 290 2'",

        # "'ba 300 2'",
        # "'ba 350 2'",
        # "'ba 400 2'",
        # "'ba 450 2'",

        # "'ba 500 2'",
        # "'ba 1000 2'",
        # "'ba 2000 2'",
        # "'ba 4000 2'",
        # "'ba 8000 2'",
        # "'ba 16000 2'",
        # "'ba 32000 2'",
        # "'ba 64000 2'",
        # "'ba 128000 2'",
        # "'ba 256000 2'",
        # "'ba 512000 2'",
        # "'ba 1024000 2'",
        # "'ba 2048000 2'",
        # "'ba 4096000 2'",
        # "'ba 8192000 2'",
        # "'ba 16384000 2'",
        # "'ba 32768000 2'",
        # "'ba 65536000 2'",
        # "'ba 131072000 2'",

        # "'flower 4 2 2'",
        # "'flower 12 2 2'",
        "'flower 44 2 2'",
        "'flower 172 2 2'",
        "'flower 684 2 2'",
        "'flower 2732 2 2'",
        "'flower 10924 2 2'",
        "'flower 43692 2 2'",
        "'flower 174764 2 2'",
        "'flower 699052 2 2'",
        "'flower 2796204 2 2'",
        # "'flower 11184812 2 2'",
    ]

    methods = [
        "sketch",
        # "coloring",
        # "memb",
        # "burning",
        # "cbb",
    ]

    exp_tag = "coverage-1.0-large-graph-scalability"
    rad_max = 10000

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
                ["qsub", "-l", "walltime=48:00:00,nodes=node06:ppn=24", "-N", job_name], stdin=p1.stdout)
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
