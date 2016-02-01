#!/usr/bin/env python
import subprocess


def torque_nageru():
    graph_names = [
        # /law.di.unimi.it/
        # "uk-2002.agl",
        # "hollywood-2011.agl",
        # "indochina-2004.agl",
        # "hollywood-2009.agl",
        # "enwiki-2013.agl",
        # "ljournal-2008.agl",
        # "uk-2007-05@1000000.agl",
        # "dewiki-2013.agl",
        # "frwiki-2013.agl",
        # "itwiki-2013.agl",
        # "eswiki-2013.agl",
        # "eu-2005.agl",
        # "in-2004.agl",
        # "dblp-2011.agl",
        # "amazon-2008.agl",
        # "cnr-2000.agl",

        # /snap.stanford.edu/
        "com-lj.ungraph.agl",
        "soc-pokec-relationships.agl",
        "cit-Patents.agl",
        "as-skitter.agl",
        "gplus_combined.agl",
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
    ]

    exp_tag = "real"
    sketch_k = 128
    pass_num = 1000
    upper_param = 1.0

    for rad in range(1, 10):
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
