#!/usr/bin/env python
import subprocess


def torque_nageru():
    graph_names = [

        # "web-Stanford.agl",
        # "web-NotreDame.agl",
        # "twitter_combined.agl",
        # # "amazon0302.agl",
        # "com-dblp.ungraph.agl",
        # "com-amazon.ungraph.agl",
        # "soc-Slashdot0902.agl",
        # # "soc-Slashdot0811.agl",
        # "email-EuAll.agl",
        # "soc-Epinions1.agl",
        # "cit-HepPh.agl",
        # "email-Enron.agl",
        # "ca-AstroPh.agl",
        # "cit-HepTh.agl",
        # "p2p-Gnutella31.agl",
        # "ca-HepPh.agl",
        # "ca-CondMat.agl",
        # # "p2p-Gnutella30.agl",
        # "as-caida20071105.agl",
        # # "p2p-Gnutella24.agl",
        # "wiki-Vote.agl",
        # # "p2p-Gnutella25.agl",
        # "facebook_combined.agl",
        # "ca-HepTh.agl",
        # # "p2p-Gnutella04.agl",
        # "oregon2_010526.agl",
        # # "oregon2_010519.agl",
        # # "oregon2_010512.agl",
        # # "oregon2_010414.agl",
        # # "oregon2_010421.agl",
        # # "oregon2_010428.agl",
        # # "oregon2_010505.agl",
        # # "oregon2_010331.agl",
        # # "oregon2_010407.agl",
        # # "p2p-Gnutella05.agl",
        # # "p2p-Gnutella06.agl",
        # # "oregon1_010526.agl",
        # # "oregon1_010519.agl",
        # # "oregon1_010512.agl",
        # # "oregon1_010505.agl",
        # # "oregon1_010421.agl",
        # # "oregon1_010428.agl",
        # # "oregon1_010414.agl",
        # # "oregon1_010407.agl",
        # # "oregon1_010331.agl",
        # # "p2p-Gnutella09.agl",
        # "ca-GrQc.agl",
        # # "p2p-Gnutella08.agl",
        # "as20000102.agl",
        # # "com-friendster.ungraph.agl",
        # "com-orkut.ungraph.agl",
        # "soc-LiveJournal1.agl",
        # "com-lj.ungraph.agl",
        # "soc-pokec-relationships.agl",
        # "cit-Patents.agl",
        # "as-skitter.agl",
        # "gplus_combined.agl",
        # "wiki-Talk.agl",
        # "roadNet-CA.agl",
        # "web-BerkStan.agl",
        # "web-Google.agl",
        # "roadNet-TX.agl",
        # "roadNet-PA.agl",
        # "com-youtube.ungraph.agl",
        # "amazon0601.agl",
        # "amazon0505.agl",
        # "amazon0312.agl",
    ]

    exp_tag = "snap"
    sketch_k = 128
    pass_num = 1000
    upper_param = 1.0
    coverage = 1.0

    for graph_name in graph_names:
        command = "/home/kenkoooo/fractal-dimension/agl/bin/box_cover --force_undirected "
        command = command + " --type agl "
        command = command + " --graph /data/snap.stanford.edu/" + graph_name

        command = command + " --final_coverage " + str(coverage)
        command = command + " --sketch_k " + str(sketch_k)
        command = command + " --pass " + str(pass_num)
        command = command + " --rad_min " + str(0)
        command = command + " --rad_max " + str(10000)
        command = command + " --upper_param " + str(upper_param)
        command = command + " --exp_tag " + exp_tag
        job_name = exp_tag + "-sketch-" + graph_name.replace("'", "").replace(" ", "-") +\
            "-k." + str(sketch_k) + \
            "-rad." + str(10000) + \
            "-upper_param." + str(upper_param)
        p1 = subprocess.Popen(["echo", command], stdout=subprocess.PIPE)
        p2 = subprocess.Popen(
            ["qsub", "-l", "walltime=48:00:00,nodes=long:ppn=12", "-p", "-10", "-N", job_name], stdin=p1.stdout)
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
