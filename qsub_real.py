#!/usr/bin/env python
import subprocess


def torque_nageru():
    graph_names = ["uk-2002.agl", "hollywood-2011.agl",
                   "enwiki-2013.agl", "ljournal-2008.agl", "arabic-2005.agl", ]

    exp_tag = "96G_real"
    sketch_k = 128
    pass_num = 1000
    upper_param = 1.0
    coverage = 1.0

    for graph_name in graph_names:
        for rad in xrange(1, 35):
            command = "/home/kenkoooo/fractal-dimension/agl/bin/box_cover --force_undirected "
            command = command + " --type agl "
            command = command + " --graph /data/law.di.unimi.it/" + graph_name

            command = command + " --final_coverage " + str(coverage)
            command = command + " --sketch_k " + str(sketch_k)
            command = command + " --pass " + str(pass_num)
            command = command + " --rad_min " + str(rad)
            command = command + " --rad_max " + str(rad)
            command = command + " --upper_param " + str(upper_param)
            command = command + " --exp_tag " + exp_tag
            job_name = exp_tag + "-sketch-" + graph_name.replace("'", "").replace(" ", "-") +\
                "-k." + str(sketch_k) + \
                "-rad." + str(rad) + \
                "-upper_param." + str(upper_param)
            p1 = subprocess.Popen(["echo", command], stdout=subprocess.PIPE)
            p2 = subprocess.Popen(
                ["qsub", "-l", "walltime=48:00:00,nodes=long:ppn=24", "-p", "-10", "-N", job_name], stdin=p1.stdout)
            p1.stdout.close()
            output = p2.communicate()[0]

    # # MEMB
    # for rad in range(1, 26):
    #     for graph_name in graph_names:
    #         command = "/home/kenkoooo/fractal-dimension/agl/bin/box_cover --force_undirected --type agl"
    #         command = command + " --graph /data/snap.stanford.edu/" + graph_name
    #         command = command + " --rad_min " + str(rad)
    #         command = command + " --rad_max " + str(rad)
    #         command = command + " --method " + "memb"
    #         command = command + " --exp_tag " + exp_tag
    #         # method_name = "memb"
    #         # command = command + " --method " + method_name
    #         # job_name = method_name + "-" + graph_name.replace("'", "").replace(" ", "-")
    #         job_name = exp_tag + "-memb-" + \
    #             graph_name.replace("'", "").replace(
    #                 " ", "-") + "-rad." + str(rad)
    #         p1 = subprocess.Popen(["echo", command], stdout=subprocess.PIPE)
    #         p2 = subprocess.Popen(
    #             ["qsub", "-l", "walltime=48:00:00", "-N", job_name], stdin=p1.stdout)
    #         p1.stdout.close()
    #         output = p2.communicate()[0]

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
