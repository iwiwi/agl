#!/usr/bin/env python
import subprocess


def torque_nageru():
    sketch_ks = [128]
    pass_nums = [1000]
    upper_params = [1.0]
    graph_names = [
        "'flower 42 3 4'",
        "'flower 287 3 4'",
        "'flower 2002 3 4'",
        "'flower 14007 3 4'",
        "'flower 98042 3 4'",
        "'flower 686287 3 4'",
        "'flower 4804002 3 4'",
        "'flower 33628007 3 4'"
    ]

    exp_tag = "scalability"
    rad_max = 1000000

    # for upper_param in upper_params:
    #     for sketch_k in sketch_ks:
    #         for pass_num in pass_nums:
    #             for graph_name in graph_names:
    #                 command = "/home/kenkoooo/fractal-dimension/agl/bin/box_cover --force_undirected "

    #                 # command = command + " --type agl "
    #                 # command = command + " --graph /data/snap.stanford.edu/" + graph_name

    #                 command = command + " --type gen "
    #                 command = command + " --graph " + graph_name
    #                 command = command + " --rad_analytical "

    #                 command = command + " --sketch_k " + str(sketch_k)
    #                 command = command + " --pass " + str(pass_num)
    #                 command = command + " --rad_max " + str(rad_max)
    #                 command = command + " --upper_param " + str(upper_param)
    #                 command = command + " --exp_tag " + exp_tag
    #                 job_name = exp_tag + "-sketch-" + graph_name.replace("'", "").replace(" ", "-") +\
    #                     "-k." + str(sketch_k) + \
    #                     "-pass." + str(pass_num) + \
    #                     "-upper_param." + str(upper_param)
    #                 p1 = subprocess.Popen(
    #                     ["echo", command], stdout=subprocess.PIPE)
    #                 p2 = subprocess.Popen(
    #                     ["qsub", "-l", "walltime=24:00:00", "-N", job_name], stdin=p1.stdout)
    #                 p1.stdout.close()
    #                 output = p2.communicate()[0]
    # # MEMB
    for graph_name in graph_names:
        command = "/home/kenkoooo/fractal-dimension/agl/bin/box_cover --force_undirected"
        command = command + " --type gen "
        command = command + " --graph " + graph_name
        command = command + " --rad_analytical "
        command = command + " --rad_max " + str(rad_max)
        command = command + " --method " + "coloring"
        command = command + " --exp_tag " + exp_tag

        job_name = exp_tag + "-coloring-" + \
            graph_name.replace("'", "").replace(" ", "-")
        p1 = subprocess.Popen(["echo", command], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(
            ["qsub", "-l", "walltime=24:00:00", "-N", job_name], stdin=p1.stdout)
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
