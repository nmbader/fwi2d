#!/usr/bin/env python3
import os
import sys
import argparse
import seppy

def parse_args():
  args=[]
  files=[]
  for arg in sys.argv[1:]:
    if arg =="-f":  args.append("force")
    elif arg =="-i":  args.append("interactive")
    elif arg =="-q":  args.append("quiet")
    elif arg =="-v":  args.append("verb")
    else: files.append(arg)
  if files == None: print("Nothing to be deleted")
  return args,files

def remove_file(file):
  sep = seppy.sep()
  faxes = sep.read_header(file)
  try:
    os.remove(sep.hdict["in"])
  except OSError:
    pass
  os.remove(file)

if __name__ == "__main__":
  
  args,files=parse_args()
  for file in files:
    if os.path.exists(file):
      remove_file(file)
  sys.exit(0)