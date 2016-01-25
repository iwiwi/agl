#!/usr/bin/env python
import subprocess


def torque_nageru():
    graph_names = [
        # 13895
        # "as20000102.agl",

        # 20777
        # "p2p-Gnutella08.agl",

        # 22000
        # "oregon1_010407.agl",

        # 22003
        # "oregon1_010331.agl",

        # 22470
        # "oregon1_010414.agl",

        # 22494
        # "oregon1_010428.agl",

        # 22608
        # "oregon1_010505.agl",

        # 22678
        # "oregon1_010512.agl",

        # 22725
        # "oregon1_010519.agl",

        # 22748
        # "oregon1_010421.agl",

        # 23410
        # "oregon1_010526.agl",

        # 26013
        # "p2p-Gnutella09.agl",

        # 28980
        # "ca-GrQc.agl",

        # 30856
        # "oregon2_010407.agl",

        # 30944
        # "oregon2_010505.agl",

        # 31181
        # "oregon2_010331.agl",

        # 31304
        # "oregon2_010512.agl",

        # 31435
        # "oregon2_010428.agl",

        # 31525
        # "p2p-Gnutella06.agl",

        # 31539
        # "oregon2_010421.agl",

        # 31762
        # "oregon2_010414.agl",

        # 31839
        # "p2p-Gnutella05.agl",

        # 32288
        # "oregon2_010519.agl",

        # 32731
        # "oregon2_010526.agl",

        # 39994
        # "p2p-Gnutella04.agl",

        # 51971
        # "ca-HepTh.agl",

        # 54705
        # "p2p-Gnutella25.agl",

        # 65369
        # "p2p-Gnutella24.agl",

        # 88234
        # "facebook_combined.agl",

        # 88328
        # "p2p-Gnutella30.agl",

        # 103689
        # "wiki-Vote.agl",

        # 106762
        # "as-caida20071105.agl",

        # 147892
        # "p2p-Gnutella31.agl",

        # 186936
        # "ca-CondMat.agl",

        # 237010
        # "ca-HepPh.agl",

        # 352807
        # "cit-HepTh.agl",

        # 367662
        # "email-Enron.agl",

        # 396160
        # "ca-AstroPh.agl",

        # 420045
        # "email-EuAll.agl",

        # 421578
        # "cit-HepPh.agl",

        # 508837
        # "soc-Epinions1.agl",

        # 905468
        # "soc-Slashdot0811.agl",

        # 925872
        # "com-amazon.ungraph.agl",

        # 948464
        # "soc-Slashdot0902.agl",

        # 1049866
        # "com-dblp.ungraph.agl",

        # 1234877
        "amazon0302.agl",

        # 1497134
        "web-NotreDame.agl",

        # 1768149
        "twitter_combined.agl",

        # 2312497
        "web-Stanford.agl",

        # 2987624
        "com-youtube.ungraph.agl",

        # 3083796
        # "roadNet-PA.agl",

        # 3200440
        # "amazon0312.agl",

        # 3356824
        # "amazon0505.agl",

        # 3387388
        # "amazon0601.agl",

        # 3843320
        # "roadNet-TX.agl",

        # 5021410
        # "wiki-Talk.agl",

        # 5105039
        # "web-Google.agl",

        # 5533214
        # "roadNet-CA.agl",

        # 7600595
        # "web-BerkStan.agl",

        # 11095298
        # "as-skitter.agl",

        # 13673453
        # "gplus_combined.agl",

        # 16518948
        # "cit-Patents.agl",

        # 30622564
        # "soc-pokec-relationships.agl",

        # 34681189
        # "com-lj.ungraph.agl",

        # 68993773
        # "soc-LiveJournal1.agl",

        # 117185083
        # "com-orkut.ungraph.agl",

        # 1806067135
        # "com-friendster.ungraph.agl",


    ]

    exp_tag = "real"
    sketch_k = 128
    pass_num = 1000
    upper_param = 1.0

    for rad in range(1, 26):
        for graph_name in graph_names:
            command = "/home/kenkoooo/fractal-dimension/agl/bin/box_cover --force_undirected "

            command = command + " --type agl "
            command = command + " --graph /data/snap.stanford.edu/" + graph_name

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
            p1 = subprocess.Popen(
                ["echo", command], stdout=subprocess.PIPE)
            p2 = subprocess.Popen(
                ["qsub", "-l", "walltime=48:00:00", "-N", job_name], stdin=p1.stdout)
            p1.stdout.close()
            output = p2.communicate()[0]
    # MEMB
    for rad in range(1, 26):
        for graph_name in graph_names:
            command = "/home/kenkoooo/fractal-dimension/agl/bin/box_cover --force_undirected --type agl"
            command = command + " --graph /data/snap.stanford.edu/" + graph_name
            command = command + " --rad_min " + str(rad)
            command = command + " --rad_max " + str(rad)
            command = command + " --method " + "memb"
            command = command + " --exp_tag " + exp_tag
            # method_name = "memb"
            # command = command + " --method " + method_name
            # job_name = method_name + "-" + graph_name.replace("'", "").replace(" ", "-")
            job_name = exp_tag + "-memb-" + \
                graph_name.replace("'", "").replace(
                    " ", "-") + "-rad." + str(rad)
            p1 = subprocess.Popen(["echo", command], stdout=subprocess.PIPE)
            p2 = subprocess.Popen(
                ["qsub", "-l", "walltime=48:00:00", "-N", job_name], stdin=p1.stdout)
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
