#!/usr/bin/env python
import subprocess


def torque_nageru():
    sketch_ks = [
        # 16,
        # 32,
        # 64,
        # 128,
        # 256,
        # 512,
        1024,
    ]
    pass_nums = [1000]
    upper_params = [
        # 0.125,
        # 0.25,
        # 0.5,
        1.0,
        # 2.0,
        # 4.0,
        # 8.0,
    ]
    graph_names = [
        # "'flower 265722 1 2'",
        # "'flower 699052 1 3'",
        # "'flower 292970 1 4'",
        # "'flower 699052 2 2'",
        # "'flower 292970 2 3'",
        # "'flower 223950 2 4'",
        # "'flower 686287 3 4'",
        # "'flower 223950 3 3'",
        "'shm 312501 5 2'",
        # "'shm 67229 5 3'",
        # "'shm 390626 6 2'",
        # "'shm 156251 3 2'",
        # "'shm 234376 4 2'",
        # "'shm 470597 5 3'",
        # "'shm 588246 6 3'",
        # "'shm 468751 7 2'",

        # "'flower 29526 1 2'",
        # "'flower 10924 1 3'",
        # "'flower 11720 1 4'",
        # "'flower 58595 1 4'",
        # "'flower 11720 2 3'",
        # "'flower 58595 2 3'",
        # "'flower 37326 2 4'",
        # "'flower 37326 3 3'",
        # "'flower 14007 3 4'",
        # "'shm 12501 5 2'",
        # "'shm 15626 6 2'",
        # "'shm 46876 4 2'",
        # "'shm 31251 3 2'",
        # "'shm 67229 5 3'",
        # "'shm 12006 6 3'",
        # "'shm 18751 7 2'",
        # "'flower 88575 1 2'",
        # "'flower 43692 2 2'",

    ]

    exp_tag = "coverage-1.0-error-by-k"
    rad_max = 10000000
    coverage = 1.0

    for upper_param in upper_params:
        for sketch_k in sketch_ks:
            for pass_num in pass_nums:
                for graph_name in graph_names:
                    for x in xrange(0, 7):
                        command = "/home/kenkoooo/fractal-dimension/agl/bin/box_cover --force_undirected "

                        command = command + " --type gen "
                        command = command + " --graph " + graph_name

                        command = command + " --rad_analytical "
                        command = command + \
                            " --final_coverage " + str(coverage)

                        command = command + " --sketch_k " + str(sketch_k)
                        command = command + " --pass " + str(pass_num)

                        command = command + " --rad_max " + str(rad_max)
                        command = command + " --random_seed " + str(x + 500)

                        command = command + \
                            " --upper_param " + str(upper_param)
                        command = command + " --exp_tag " + exp_tag
                        job_name = exp_tag + "-sketch-" + graph_name.replace("'", "").replace(" ", "-") +\
                            "-k." + str(sketch_k) + \
                            "-pass." + str(pass_num) + \
                            "-coverage." + str(coverage) + \
                            "-upper_param." + str(upper_param)
                        p1 = subprocess.Popen(
                            ["echo", command], stdout=subprocess.PIPE)
                        p2 = subprocess.Popen(
                            ["qsub", "-l", "walltime=24:00:00,nodes=1:ppn=8", "-N", job_name], stdin=p1.stdout)
                        p1.stdout.close()
                        output = p2.communicate()[0]
    # MEMB
    # for graph_name in graph_names:
    #     command = "/home/kenkoooo/fractal-dimension/agl/bin/box_cover --force_undirected --type agl"
    #     command = command + " --graph /data/snap.stanford.edu/" + graph_name
    #     command = command + " --rad_max " + str(10000)
    #     command = command + " --method " + "memb"
    #     # method_name = "memb"
    #     # command = command + " --method " + method_name
    #     # job_name = method_name + "-" + graph_name.replace("'", "").replace(" ", "-")
    #     job_name = "memb-" + graph_name.replace("'", "").replace(" ", "-")
    #     p1 = subprocess.Popen(["echo", command], stdout=subprocess.PIPE)
    #     p2 = subprocess.Popen(
    #         ["qsub", "-l", "walltime=24:00:00,nodes=1:ppn=12", "-N", job_name], stdin=p1.stdout)
    #     p1.stdout.close()
    #     output = p2.communicate()[0]

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