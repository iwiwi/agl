#!/usr/bin/env python
import subprocess


def torque_nageru():
    sketch_ks = [
        128,
        # 32, 64, 256, 512
    ]
    pass_nums = [1000]
    upper_params = [
        1.0,
        0.25, 0.5,
        2.0,
        4.0,
    ]
    graph_names = [
        # "com-friendster.ungraph.agl",
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
        # "web-Stanford.agl",
        # "web-NotreDame.agl",
        # "twitter_combined.agl",
        # "amazon0302.agl",
        # "com-dblp.ungraph.agl",
        # "com-amazon.ungraph.agl",
        # "soc-Slashdot0902.agl",
        # "soc-Slashdot0811.agl",
        # "email-EuAll.agl",
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

        "'flower 265722 1 2'",
        "'flower 699052 1 3'",
        "'flower 292970 1 4'",
        "'flower 699052 2 2'",
        "'flower 292970 2 3'",
        "'flower 223950 2 4'",
        "'flower 223950 3 3'",
        "'flower 686287 3 4'",
        "'shm 312501 5 2'",
        "'shm 390626 6 2'"
    ]
    for upper_param in upper_params:
        for sketch_k in sketch_ks:
            for pass_num in pass_nums:
                for graph_name in graph_names:
                    command = "/home/kenkoooo/fractal-dimension/agl/bin/box_cover --force_undirected "

                    # command = command + " --type agl "
                    # command = command + " --graph /data/snap.stanford.edu/" + graph_name

                    command = command + " --type gen "
                    command = command + " --graph " + graph_name
                    command = command + " --rad_analytical "

                    command = command + " --sketch_k " + str(sketch_k)
                    command = command + " --pass " + str(pass_num)
                    command = command + " --rad_max " + str(1000000)
                    command = command + " --upper_param " + str(upper_param)
                    job_name = "exp01-sketch-" + graph_name.replace("'", "").replace(" ", "-") +\
                        "-k." + str(sketch_k) + \
                        "-pass." + str(pass_num) + \
                        "-upper_param." + str(upper_param)
                    p1 = subprocess.Popen(
                        ["echo", command], stdout=subprocess.PIPE)
                    p2 = subprocess.Popen(
                        ["qsub", "-l", "walltime=24:00:00", "-N", job_name], stdin=p1.stdout)
                    p1.stdout.close()
                    output = p2.communicate()[0]
    # MEMB
    # for rad_max in rad_maxs:
    #     for graph_name in graph_names:
    #         command = "/home/kenkoooo/fractal-dimension/agl/bin/box_cover --force_undirected --type gen"
    #         command = command + " --graph " + graph_name
    #         command = command + " --rad_max " + str(rad_max)
    #         command = command + " --method " + "memb"
    #         command = command + " --rad_analytical "
    #         # method_name = "memb"
    #         # command = command + " --method " + method_name
    #         # job_name = method_name + "-" + graph_name.replace("'", "").replace(" ", "-")
    #         job_name = "memb-" + graph_name.replace("'", "").replace(" ", "-")
    #         p1 = subprocess.Popen(["echo", command], stdout=subprocess.PIPE)
    #         p2 = subprocess.Popen(
    #             ["qsub", "-l", "walltime=24:00:00,nodes=1:ppn=4", "-N", job_name], stdin=p1.stdout)
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
