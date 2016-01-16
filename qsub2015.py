#!/usr/bin/env python
import subprocess

if __name__ == "__main__":
    sketch_ks = [32, 64, 128, 256, 512]
    rad_maxs = [512]
    pass_nums = [1000]
    size_uppers = [2500000, 5000000, 10000000, 20000000, 40000000]
    graph_names = [
        "'flower 88575 1 2'", "'flower 43692 1 3'", "'flower 58595 1 4'", "'flower 37326 1 5'",
        "'flower 265722 1 2'", "'flower 699052 1 3'", "'flower 292970 1 4'", "'flower 223950 1 5'",
        "'flower 699052 2 2'", "'flower 292970 2 3'", "'flower 223950 2 4'",
        "'flower 223950 3 3'", "'flower 686287 3 4'",
        "'shm 312501 5 2'",  "'shm 390626 6 2'"
    ]
    for size_upper in size_uppers:
        for sketch_k in sketch_ks:
            for rad_max in rad_maxs:
                for pass_num in pass_nums:
                    for graph_name in graph_names:

                        command = "/home/kenkoooo/fractal-dimension/agl/bin/box_cover --force_undirected --type gen"
                        command = command + " --graph " + graph_name
                        command = command + " --rad_max " + str(rad_max)
                        command = command + " --sketch_k " + str(sketch_k)
                        command = command + " --pass " + str(pass_num)
                        command = command + " --rad_analytical "
                        command = command + " --size_upper " + str(size_upper)
                        # method_name = "memb"
                        # command = command + " --method " + method_name
                        # job_name = method_name + "-" + graph_name.replace("'", "").replace(" ", "-")
                        job_name = "sketch-" + graph_name.replace("'", "").replace(" ", "-") + \
                            "-k." + str(sketch_k) + "-rad." + \
                            str(rad_max) + "-pass." + str(pass_num) + \
                            "-size_upper." + str(size_upper)
                        p1 = subprocess.Popen(
                            ["echo", command], stdout=subprocess.PIPE)
                        p2 = subprocess.Popen(
                            ["qsub", "-l", "walltime=48:00:00", "-N", job_name], stdin=p1.stdout)
                        p1.stdout.close()
                        output = p2.communicate()[0]
    # MEMB
    for rad_max in rad_maxs:
        for graph_name in graph_names:
            command = "/home/kenkoooo/fractal-dimension/agl/bin/box_cover --force_undirected --type gen"
            command = command + " --graph " + graph_name
            command = command + " --rad_max " + str(rad_max)
            command = command + " --method " + "memb"
            command = command + " --rad_analytical "
            # method_name = "memb"
            # command = command + " --method " + method_name
            # job_name = method_name + "-" + graph_name.replace("'", "").replace(" ", "-")
            job_name = "memb-" + graph_name.replace("'", "").replace(" ", "-")
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
