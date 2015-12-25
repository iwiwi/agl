import subprocess

if __name__ == "__main__":
    sketch_ks = [128, 256, 512, 1024]
    rad_maxs = [8, 16, 32, 64, 128, 256]
    pass_nums = [1, 100]
    graph_names = ["'flower 1000 2 2'", "'flower 10000 2 2'", "'flower 1000 1 2'",
                   "'flower 10000 1 2'", "'shm 1000 5 2'", "'shm 10000 5 2'"]
    for sketch_k in sketch_ks:
        for rad_max in rad_maxs:
            for pass_num in pass_nums:
                for graph_name in graph_names:

                    command = "/home/kenkoooo/fractal-dimension/agl/bin/box_cover --force_undirected --type gen"
                    command = command + " --graph " + graph_name
                    command = command + " --rad_max " + str(rad_max)
                    command = command + " --sketch_k " + str(sketch_k)
                    command = command + " --pass " + str(pass_num)
                    # method_name = "memb"
                    # command = command + " --method " + method_name
                    # job_name = method_name + "-" + graph_name.replace("'", "").replace(" ", "-")
                    job_name = "sketch-" + "k." + str(sketch_k) + "-rad." + str(rad_max) + \
                        "-pass." + str(pass_num) + graph_name.replace("'", "").replace(" ", "-")
                    p1 = subprocess.Popen(["echo", command], stdout=subprocess.PIPE)
                    p2 = subprocess.Popen(["qsub", "-l", "walltime=240:00:00", "-N", job_name], stdin=p1.stdout)
                    p1.stdout.close()
                    output = p2.communicate()[0]

    for rad_max in rad_maxs:
        for graph_name in graph_names:
            command = "/home/kenkoooo/fractal-dimension/agl/bin/box_cover --force_undirected --type gen"
            command = command + " --graph " + graph_name
            command = command + " --rad_max " + str(rad_max)
            command = command + " --method " + "memb"
            # method_name = "memb"
            # command = command + " --method " + method_name
            # job_name = method_name + "-" + graph_name.replace("'", "").replace(" ", "-")
            job_name = "memb-" + "k." + str(sketch_k) + "-rad." + str(rad_max) + \
                "-pass." + str(pass_num) + graph_name.replace("'", "").replace(" ", "-")
            p1 = subprocess.Popen(["echo", command], stdout=subprocess.PIPE)
            p2 = subprocess.Popen(["qsub", "-l", "walltime=240:00:00", "-N", job_name], stdin=p1.stdout)
            p1.stdout.close()
            output = p2.communicate()[0]
