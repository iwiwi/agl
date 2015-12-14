import sys,os,subprocess,glob,re

directory = "exec_test"

for name in ["dolphin", "karate_club"]:
  command = "bin/mincut_query -type=built_in -graph=%s --jlog_out=%s" % (name, directory)

  p = subprocess.Popen(command.split(), stdin=subprocess.PIPE, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
  p.wait()
  err = p.stderr.readlines()
  g = re.search("JLOG: (.*)", err[0])
  filename = g.group(1)
  print g
  print filename
  os.rename(directory + "/" + filename,directory + "/" + name + ".jlog")
