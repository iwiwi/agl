#!/usr/bin/env python
import subprocess


def torque_nageru():
    sketch_k = 128
    pass_num = 1000
    upper_param = 1.0
    graph_names = [

        "twitter_combined.agl",
        "amazon0302.agl",
        "com-dblp.ungraph.agl",
        # "soc-Epinions1.agl",
        # "cit-HepPh.agl",
        # "email-Enron.agl",
        # "ca-AstroPh.agl",
        # "cit-HepTh.agl",
        # "p2p-Gnutella31.agl",
        # "ca-HepPh.agl",
        # "ca-CondMat.agl",
        # "p2p-Gnutella30.agl",
        # "as-caida20071105.agl",
        # "p2p-Gnutella24.agl",
        # "wiki-Vote.agl",
        # "p2p-Gnutella25.agl",
        # "facebook_combined.agl",
        # "ca-HepTh.agl",
        # "p2p-Gnutella04.agl",
        # "oregon2_010526.agl",
        # "oregon2_010519.agl",
        # "oregon2_010512.agl",
        # "oregon2_010414.agl",
        # "oregon2_010421.agl",
        # "oregon2_010428.agl",
        # "oregon2_010505.agl",
        # "oregon2_010331.agl",
        # "oregon2_010407.agl",
        # "p2p-Gnutella05.agl",
        # "p2p-Gnutella06.agl",
        # "oregon1_010526.agl",
        # "oregon1_010519.agl",
        # "oregon1_010512.agl",
        # "oregon1_010505.agl",
        # "oregon1_010421.agl",
        # "oregon1_010428.agl",
        # "oregon1_010414.agl",
        # "oregon1_010407.agl",
        # "oregon1_010331.agl",
        # "p2p-Gnutella09.agl",
        # "ca-GrQc.agl",
        # "p2p-Gnutella08.agl",
        # "as20000102.agl",
        # "web-BerkStan.agl",
        # "web-Google.agl",
        # # "web-NotreDame.agl",
        # "web-Stanford.agl",
    ]

    methods = [
        # "sketch",
        "coloring",
        # "memb",
        "cbb",
        "burning",
    ]

    exp_tag = "main_table-real-coverage.1.0"
    rad_max = 20

    for method in methods:
        for graph_name in graph_names:
            # for rad in range(1, 10):  # 1<=rad<10
            command = "/home/kenkoooo/fractal-dimension/agl/bin/box_cover --force_undirected"
            command = command + " --type agl "
            command = command + " --graph /data/snap.stanford.edu/" + graph_name
            command = command + " --final_coverage " + str(1.0)
            command = command + " --rad_min " + str(1)
            command = command + " --rad_max " + str(rad_max)
            command = command + " --method " + method
            command = command + " --exp_tag " + exp_tag
            command = command + " --sketch_k " + str(sketch_k)
            command = command + " --pass " + str(pass_num)
            command = command + " --upper_param " + str(upper_param)

            job_name = exp_tag + "-" + method + "-" + \
                graph_name.replace("'", "").replace(" ", "-") +\
                "-rad." + str(rad_max)
            p1 = subprocess.Popen(
                ["echo", command], stdout=subprocess.PIPE)
            p2 = subprocess.Popen(
                ["qsub", "-l", "walltime=24:00:00", "-N", job_name], stdin=p1.stdout)
            p1.stdout.close()
            output = p2.communicate()[0]

if __name__ == "__main__":
    torque_nageru()
