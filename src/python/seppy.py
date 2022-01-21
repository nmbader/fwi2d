# Author Joseph Jennings

import numpy as np
import os as opsys
import socket
import getpass
import __main__ as main
import datetime
import string
import random


class axes:
  """ Axes of regularly sampled data"""

  def __init__(self, n, o, d, label=None):
    self.ndims = len(n)
    if len(o) != self.ndims or len(d) != self.ndims:
      raise ValueError("the size of n, o and d do not match.")
    self.n = n
    self.o = o
    self.d = d
    self.label = label

  def get_nelem(self):
    return np.prod(self.n)

  def add_axis(self, nin, oin, din):
    self.n.append(nin)
    self.o.append(oin)
    self.d.append(din)
    self.ndims += 1


class sep:
  """ Utility for reading, writing and plotting SEP files """

  def __init__(self, argv=None):
    self.hdict = {}
    self.haxes = {}
    self.argv = argv
    self.hostname = self.gethostname()

  def get_fname(self, tag):
    """ Gets the file name with the associated tag """
    if (self.argv is None):
      raise Exception("Must pass sys.argv to constructor in order to use tags")
    fname = None
    for iarg in self.argv:
      keyval = iarg.split('=')
      if (len(keyval) != 2):
        continue
      if (keyval[0] == tag):
        fname = keyval[1]
        break

    return fname

  def read_header(self, ifname, tag=None, hdict=False):
    """ Reads a SEP header file from tag and returns the axes """
    self.hdict = {}
    # Get the filename with the given tag
    if tag is None and ifname is None:
      raise ValueError("Need either a tag or inputfile for reading")
    if (tag is not None):
      hin = self.get_fname(tag)
      assert (hin is not None), "No file associated with tag '%s'" % (tag)
    else:
      hin = ifname
    # Work with tags if argv has been passed
    if self.argv is not None:
      hout = self.get_fname("out")
      fout = None
      if hout is not None:
        fout = open(hout, "w")
    # First read header into dictionary
    for line in open(hin).readlines():
      # If output file, copy history to output
      if tag == "in" and hout is not None:
        fout.write(line)
      splitspace = line.split(' ')
      for item in splitspace:
        spliteq = item.split('=')
        if (len(spliteq) == 1):
          continue
        spliteq[0] = spliteq[0].replace('\n', '')
        spliteq[0] = spliteq[0].replace('\t', '')
        spliteq[1] = spliteq[1].replace('\n', '')
        spliteq[1] = spliteq[1].replace('\t', '')
        self.hdict[spliteq[0]] = spliteq[1]
    # Check if it found a binary
    assert (
        "in" in self.hdict
    ), "Error: header in file %s does not have an associated binary." % (hin)
    self.hdict["in"] = self.hdict["in"].replace('"', '')
    # Read the header info into a list of axes
    ns = []
    os = []
    ds = []
    lbls = []
    for n in range(1, 7):
      nkey = "n" + str(n)
      okey = "o" + str(n)
      dkey = "d" + str(n)
      lblkey = "label" + str(n)
      if n == 1:
        assert (nkey
                in self.hdict), "Error: header in file %s has no n1." % (hin)
      if nkey in self.hdict and okey in self.hdict and dkey in self.hdict:
        ns.append(int(self.hdict[nkey]))
        os.append(float(self.hdict[okey]))
        ds.append(float(self.hdict[dkey]))
        if (lblkey in self.hdict):
          lbls.append(self.hdict[lblkey])
        else:
          lbls.append(" ")

    # Remove ones at the end
    for n in ns:
      if (ns[-1] == 1.0 and len(ns) != 1):
        del ns[-1]
        del os[-1]
        del ds[-1]
        del lbls[-1]

    # Take care of the remaining
    if (ns[-1] == 1.0 and len(ns) != 1):
      del ns[-1]
      del os[-1]
      del ds[-1]
      del lbls[-1]

    if (lbls == []):
      if (hdict):
        return self.hdict, axes(ns, os, ds, None)
      else:
        return axes(ns, os, ds, None)
    else:
      if (hdict):
        return self.hdict, axes(ns, os, ds, lbls)
      else:
        return axes(ns, os, ds, lbls)

  def read_header_dict(self, ifname, tag=None):
    """ Reads a SEP header file and returns a dictionary """
    hdict = {}
    # Get the filename with the given tag
    if tag is None and ifname is None:
      raise ValueError("Need either a tag or inputfile for reading")
    if tag is not None:
      hin = self.get_fname(tag=tag)
      assert (hin is not None), "No file associated with tag '%s'" % (tag)
    else:
      hin = ifname
    # Read header into dictionary
    for line in open(hin).readlines():
      splitspace = line.split(' ')
      for item in splitspace:
        spliteq = item.split('=')
        if (len(spliteq) == 1):
          continue
        spliteq[0] = spliteq[0].replace('\n', '')
        spliteq[0] = spliteq[0].replace('\t', '')
        spliteq[1] = spliteq[1].replace('\n', '')
        spliteq[1] = spliteq[1].replace('\t', '')
        hdict[spliteq[0]] = spliteq[1]

    return hdict

  def read_file(self, ifname, form='xdr', safe=False, tag=None):
    """ Reads a SEP file from tag and returns the data and the axes """
    faxes = self.read_header(ifname, tag)
    # Get the correct data type
    esize = int(self.hdict['esize'])
    dtype = self.get_dtype(form, esize)
    # Read in the file
    if (safe):
      dat = np.zeros(faxes.get_nelem())
      with open(self.hdict["in"], 'rb') as f:
        dat[:] = np.fromfile(f, dtype=dtype)
    else:
      with open(self.hdict["in"], 'rb') as f:
        dat = np.fromfile(f, dtype=dtype)

    return faxes, dat

  def read_wind(self, ifname, fw, nw, form='xdr', safe=False, tag=None):
    """
    Reads a portion of a SEP file (from last axis) and
    returns the data and the axes
    """
    faxes = self.read_header(ifname, tag)
    # Get the correct data type
    esize = int(self.hdict['esize'])
    dtype = self.get_dtype(form, esize)
    # Compute the offset and update output axes
    offset = int(esize * fw * np.prod(faxes.n[:-1]))
    faxes.n[-1] = nw
    count = faxes.get_nelem()
    if (nw == 1):
      del faxes.n[-1]
    if (safe):
      dat = np.zeros(faxes.n)
      with open(self.hdict["in"], 'rb') as f:
        dat[:] = np.fromfile(f, dtype=dtype, count=count, offset=offset)
    else:
      with open(self.hdict["in"], 'rb') as f:
        dat = np.fromfile(f, dtype=dtype, count=count, offset=offset)

    return faxes, dat

  def from_header(self, ifname, keys, tag=None):
    """ Given a list of keys (strings), returns the values from the header """
    odict = {}
    # Read the header dictionary
    thdict = self.read_header_dict(ifname, tag)
    # Loop over all keys
    for ikey in keys:
      if ikey in thdict:
        odict[ikey] = thdict[ikey]

    return odict

  def write_header(
      self,
      ofname,
      ns,
      esize,
      os=None,
      ds=None,
      ofaxes=None,
      tag=None,
      dpath=None,
      form='xdr',
  ):
    """
    Writes header information to SEP file and returns
    the path to the output
    """
    fout = None
    assert (tag is not None or ofname
            is not None), "Need a tag or output file name to write a header."
    if tag is not None:
      ofname = self.get_fname(tag)

    if (tag == "out"):
      if ofname is None:
        raise ValueError("No output file name found. \
        To work with tags you must pass argv to the constructor.")
      fout = open(ofname, "a")
    else:
      fout = open(ofname, "w+")
    # Write the first line
    fout.write('\n' + self.get_fline() + '\n')
    # Get the datapath
    if len(ofname.split('/')) > 1:
      ofname = ofname.split('/')[-1]
    opath = None
    if dpath is None:
      opath = self.get_datapath() + ofname + "@"
    else:
      opath = dpath + ofname + "@"
    fout.write('\t\tsets next: in="%s"\n' % (opath))
    if (ofaxes is not None):
      # Print axes
      for k in range(ofaxes.ndims):
        if ofaxes.label is not None:
          fout.write("\t\tn%d=%d o%d=%f d%d=%.12f label%d=%s\n" %
                     (k + 1, ofaxes.n[k], k + 1, ofaxes.o[k], k + 1,
                      ofaxes.d[k], k + 1, ofaxes.label[k]))
        else:
          fout.write(
              "\t\tn%d=%d o%d=%f d%d=%.12f\n" %
              (k + 1, ofaxes.n[k], k + 1, ofaxes.o[k], k + 1, ofaxes.d[k]))
    else:
      # Write with as much info as they provide
      ndim = len(ns)
      if (os is None):
        os = np.zeros(ndim)
      if (ds is None):
        ds = np.ones(ndim)
      else:
        if len(ds) != ndim:
          ds = np.pad(
              ds,
              (0, ndim - len(ds)),
              mode='constant',
              constant_values=1,
          )
        if len(os) != ndim:
          os = np.pad(
              os,
              (0, ndim - len(os)),
              mode='constant',
              constant_values=0.0,
          )
      for k in range(ndim):
        fout.write("\t\tn%d=%d o%d=%f d%d=%.12f\n" %
                   (k + 1, ns[k], k + 1, os[k], k + 1, ds[k]))

    if esize == 1:
      fout.write('\t\tdata_format="native_byte" esize=%d\n' % (esize))
    else:
      if (form == 'xdr'):
        fout.write('\t\tdata_format="xdr_float" esize=%d\n' % (esize))
      elif (form == 'native'):
        fout.write('\t\tdata_format="native_float" esize=%d\n' % (esize))
      else:
        print("Error: format %s not recognized" % (form))

    fout.close()

    return opath

  def get_dtype(self, form, esize):
    """ Returns the datatype based on the esize and output file format """
    # Real
    if (esize == 4):
      if (form == 'xdr'):
        dtype = '>f'
      elif (form == 'native'):
        dtype = '<f'
      else:
        raise Exception("Failed to read in file. Format %s not recognized\n" %
                        (form))
    # Complex
    elif (esize == 8):
      if (form == 'xdr'):
        dtype = '>c8'
      elif (form == 'native'):
        dtype = '<c8'
      else:
        raise Exception("Failed to read in file. Format %s not recognized\n" %
                        (form))
    elif (esize == 1):
      dtype = np.ubyte
    else:
      raise Exception(
          "Failed to read in file. Must have esize 4 or 8 (esize=%d)" % (esize))

    return dtype

  def get_esize_dtype(self, data, form):
    """ Returns the esize, type and endianness for writing to a file """
    if ("f" in "%s" % (data.dtype)):
      esize = 4
      if (form == 'xdr'):
        dtype = '>f'
      elif (form == 'native'):
        dtype = '<f'
      else:
        raise ValueError('Failed to write file. Format %s not recognized\n' %
                         (form))
    elif ("c" in "%s" % (data.dtype)):
      esize = 8
      if (form == 'xdr'):
        dtype = '>c8'
      elif (form == 'native'):
        dtype = '<c8'
      else:
        raise ValueError('Failed to write file. Format %s not recognized\n' %
                         (form))
    else:
      raise ValueError("Error: can only write real or complex data")

    return esize, dtype

  def write_file(
      self,
      ofname,
      data,
      os=None,
      ds=None,
      ofaxes=None,
      tag=None,
      dpath=None,
      form='xdr',
  ):
    """ Writes data and axes to a SEP header and binary """
    # Get data type and esize (real or complex)
    if ("f" in "%s" % (data.dtype)):
      esize = 4
      if (form == 'xdr'):
        dtype = '>f'
      elif (form == 'native'):
        dtype = '<f'
      else:
        raise ValueError('Failed to write file. Format %s not recognized\n' %
                         (form))
    elif ("c" in "%s" % (data.dtype)):
      esize = 8
      if (form == 'xdr'):
        dtype = '>c8'
      elif (form == 'native'):
        dtype = '<c8'
      else:
        raise ValueError('Failed to write file. Format %s not recognized\n' %
                         (form))

    elif (data.dtype == 'uint8'):
      esize = 1
      dtype = np.ubyte

    else:
      raise ValueError("Error: can only write real or complex data")
    # Write header
    opath = self.write_header(ofname, data.shape, esize, os, ds, ofaxes, tag,
                              dpath, form)
    with open(opath, 'wb') as f:
      data.flatten('F').astype(dtype).tofile(f)

  def append_file(
      self,
      ofname,
      data,
      os=None,
      ds=None,
      newaxis=False,
      ofaxes=None,
      tag=None,
  ):
    """
    Appends the input data to a SEP file
    (use this after first writing the file)
    """
    # Read in the header to be appended and get the last axis
    hdict, faxes = self.read_header(ofname, tag, hdict=True)
    diffn = data.ndim - faxes.ndims
    if diffn == -1:
      # Appending one new example to an already appended file
      # (output is larger than input)
      odim = faxes.ndims
      nnew = faxes.n[-1] + 1
      onew = faxes.o[-1]
      dnew = faxes.d[-1]

    elif (diffn == 0):
      if (newaxis is False):
        # Appending examples to an already appended file
        # (input and output are same size)
        odim = faxes.ndims
        nnew = faxes.n[-1] + data.shape[-1]
        onew = faxes.o[-1]
        dnew = faxes.d[-1]

      else:
        # Appending one example for the first time
        # (input and output are same size)
        odim = faxes.ndims + 1
        nnew = 2
        if (os is not None):
          if (len(os) > faxes.ndims):
            onew = os[-1]
          else:
            onew = 0.0
        else:
          onew = 0.0
        if (ds is not None):
          if (len(ds) > faxes.ndims):
            dnew = ds[-1]
          else:
            dnew = 1.0
        else:
          dnew = 1.0

    elif (diffn == 1):
      # Appending >= 1 examples to a non-appended file
      # (input is larger than output)
      odim = faxes.ndims + 1
      nnew = 1 + data.shape[-1]
      if (os is not None):
        onew = os[-1]
      else:
        onew = 0.0
      if (ds is not None):
        dnew = ds[-1]
      else:
        dnew = 1.0

    else:
      raise ValueError("Input shape not correct for appending to file. \
                       Output ndim is %d and input is %d" %
                       (len(data.shape), len(faxes.n)))

    # Write info to header
    if ofname is None:
      ofname = self.get_fname(tag)
    fout = open(ofname, "a")
    fout.write("\n\t\tn%d=%d o%d=%f d%d=%f\n" %
               (odim, nnew, odim, onew, odim, dnew))
    fout.close()

    # Append data to binary
    if (hdict['data_format'] == '"xdr_float"'):
      form = 'xdr'
    elif (hdict['data_format'] == '"native_float"'):
      form = 'native'
    esize, dtype = self.get_esize_dtype(data, form)
    with open(hdict['in'], 'ab') as f:
      data.flatten('F').astype(dtype).tofile(f)

  def to_header(self, ofname, info, tag=None):
    """ Writes any auxiliary information to header """
    fout = None
    if tag is not None and ofname is not None:
      raise ValueError("Need a tag or output file name to write a header.")
    if tag is not None:
      ofname = self.get_fname(tag)

    # Open file header
    if tag == "out":
      if ofname is None:
        raise ValueError("No output file name found. \
                         To work with tags, you must pass \
                         argv to the constructor")
      fout = open(ofname, "a")
    else:
      fout = open(ofname, "a+")
    fout.write('\n' + info)
    fout.close()

  def write_dummyaxis(self, ofname, dim, tag=None):
    """ Writes a single axis to a SEP header """
    if ofname is None:
      ofname = self.get_fname(tag)

    fout = open(ofname, "a")
    fout.write("\n\t\tn%d=1 o%d=0.0 d%d=1.0\n" % (dim, dim, dim))
    fout.close()

  def append_to_movie(
      self,
      ofname,
      ofaxes,
      data,
      niter,
      tag=None,
      dpath=None,
      form='xdr',
  ):
    """ Appends to a file for an inversion movie"""
    if ofname is None:
      ofname = self.get_fname(tag)

    # Write an extra line to the header
    fout = open(ofname, "a")
    odim = ofaxes.ndims + 1
    fout.write("\n\t\tn%d=%d o%d=0.0 d%d=1.0\n" % (odim, niter, odim, odim))
    fout.close()

    # Append the data to the binary
    opath = None
    if dpath is None:
      opath = opsys.path.join(self.get_datapath(), ofname) + '@'
    else:
      opath = opsys.path.join(dpath, ofname) + '@'
    with open(opath, 'ab') as f:
      if (form == 'xdr'):
        data.flatten('F').astype('>f').tofile(f)
      elif (form == 'native'):
        data.flatten('F').astype('<f').tofile(f)

  def get_fline(self):
    """ Returns the first line of the program header """
    if (self.isnotebook()):
      fline = 'ipython'
    elif (self.argv is None):
      fline = "%s" % (main.__file__)
    else:
      fline = self.argv[0]
    # Get user and hostname
    username = getpass.getuser()
    fline += ":\t" + username + "@" + self.hostname + "\t\t"
    # Get time and date
    time = datetime.datetime.today()
    fline += time.strftime("%a %b %d %H:%M:%S %Y")

    return fline

  def get_datapath(self):
    """ Gets the set datpath for writing SEP binaries """
    dpath = None
    # Check if passed as argument
    if (self.argv is not None):
      dpath = self.get_fname("datapath")
    if (dpath is None):
      # Look in home directory
      datstring = opsys.environ['HOME'] + "/.datapath"
      if (opsys.path.exists(datstring) is True):
        nohost = ''
        # Assumes structure as host datapath=path
        for line in open(datstring).readlines():
          hostpath = line.split()
          if (len(hostpath) == 1):
            nohost = hostpath
          elif (len(hostpath) == 2):
            # Check if hostname matches
            if (self.hostname == hostpath[0]):
              dpath = hostpath[1].split('=')[1]
              break
        if (dpath is None and nohost is not None):
          dpath = nohost[0].split('=')[1]
    # Lastly, look at environment variable
    elif (dpath is None and "DATAPATH" in opsys.environ):
      dpath = opsys.environ['DATAPATH']
    # else:
    #   dpath = '/tmp/'

    return dpath

  def gethostname(self, alias=True):
    """
    An extension of socket.gethostname() where the hostname alias
    is returned if requested. This makes this code more
    compatible with the IO in SEPlib
    """
    hname = socket.gethostname()
    if (alias):
      with open('/etc/hosts', 'r') as f:
        for line in f.readlines():
          cols = line.lower().split()
          if (len(cols) > 1):
            if (hname == cols[1] and len(cols) > 2):
              hname = cols[2]

    return hname

  def isnotebook(self):
    """ Checks if running in a notebook/ipython """
    try:
      shell = get_ipython().__class__.__name__
      if (shell == 'ZMQInteractiveShell' or
          shell == 'TerminalInteractiveShell'):
        return True  # Jupyter notebook or IPython
      else:
        return False  # Other type (?)
    except NameError:
      return False

  def id_generator(self, size=6, chars=string.ascii_uppercase + string.digits):
    """ Creates a random string with uppercase letters and integers """
    return ''.join(random.choice(chars) for _ in range(size))

  def yn2zoo(self, yn):
    """ Converts a 'y' or 'n' to an integer """
    if (yn == "n"):
      zoo = 0
    else:
      zoo = 1

    return zoo

  def read_list(self, arg, default, dtype='int'):
    """
    Reads in comma delimited string at the command line
    into a python list
    """
    if (len(arg) == 0):
      return default
    else:
      if (dtype == 'int'):
        olist = [int(val) for val in arg.split(',')]
        return olist
      elif (dtype == 'float'):
        olist = [float(val) for val in arg.split(',')]
        return olist
      else:
        print("Type %s not recognized. Returning default")
        return default


def bytes2float(img):
  """ Converts an array of bytes (uint8) to floats """
  imgtmp = img - 255 / 2.
  imgtmp *= 1. / 255
  return imgtmp
