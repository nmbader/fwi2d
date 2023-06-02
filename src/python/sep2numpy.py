#!/usr/bin/env python3
import numpy as np
import argparse
import seppy


if __name__ == "__main__":

  parser = argparse.ArgumentParser()
  parser.add_argument("--input", type=str, default="None", help="Input file in numpy or SEPlib format")
  parser.add_argument("--output", type=str, default="None", help="Output file in numpy or SEPlib format")
  parser.add_argument("--mode", type=int, default=0, help="0: convert sep file to numpy. 1: convert numpy file to sep format")
  parser.add_argument("--datapath", type=str, default=None, help="Datapath to save binary data in mode 1")
  args = parser.parse_args()

  sep = seppy.sep()

  if (args.mode==0):
    axes, data = sep.read_file(args.input)
    npdat = data.reshape(axes.n,order='F').T
    np.save(args.output, npdat)

  elif (args.mode==1):
    npdat = np.load(args.input)
    sep.write_file(args.output, np.transpose(npdat), ds=np.ones(npdat.ndim), os=np.zeros(npdat.ndim), dpath=args.datapath)