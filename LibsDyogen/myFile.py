# -*- coding: utf-8 -*-

# LibsDyogen_py3 version 0 (2018/11/16)
# Copyright © 2015 IBENS/Dyogen : Matthieu MUFFATO, Joseph LUCAS and Hugues ROEST CROLLIUS
# Licenses GPL v3 and CeCILL v2

"""file management functions"""

import itertools
import collections
import os
import sys

null = open(os.devnull, 'w')


def codec_hook(f):
    """For Python 3: if the file is in a compressed/binary format, convert
    input from `str` to `bytes`, and vice-versa for the output."""
    if sys.version_info.major >= 3:
        orig_read = f.read
        orig_readlines = f.readlines
        orig_readline = f.readline
        orig_write = f.write
        orig_writelines = f.writelines
        def readdecode(*args, **kwargs):
            return orig_read().decode()
        def readlinesdecode(*args, **kwargs):
            return [line.decode() for line in orig_readlines()]
        def readlinedecode(*args, **kwargs):
            return orig_readline().decode()
        def writeencode(text):
            return orig_write(text.encode())
        def writelinesencode(lines):
            return orig_writelines((text.encode() for text in lines))
        f.read       = readdecode
        f.readlines  = readlinesdecode
        f.readline   = readlinedecode
        f.write      = writeencode
        f.writelines = writelinesencode
    return f


class myTSV:
    """tabular file management"""

    csvProxy = collections.namedtuple("csvProxy", ['file','csvobject'])

    @staticmethod
    def reader(fileName, **keywords):
        """read with the csv module"""
        import csv
        f = openFile(fileName, 'r')
        return myTSV.csvProxy(f,csv.reader(f, delimiter="\t", quoting=csv.QUOTE_NONE, lineterminator="\n", **keywords))

    @staticmethod
    def writer(fileName):
        """write with the csv module"""
        import csv
        f = openFile(fileName, 'w')
        return myTSV.csvProxy(f,csv.writer(f, delimiter="\t", quoting=csv.QUOTE_NONE, lineterminator="\n"))

    @staticmethod
    def printLine(line, delim = "\t", func = str):
        """return the prepared line for printing"""
        return delim.join(func(x) for x in line)


    @staticmethod
    def readTabular(filename, type_list, delim = '\t'):
        """read a tabular file, convert columns separated by delim depending on
        the type_list"""
        f = openFile(filename, 'r')
        # list of each column type
        new_type_list = []
        for x in type_list:
            if type(x) == type:
                new_type_list.append(x)
            else:
                new_type_list.extend([x[0]] * x[1])
        # run through the file (parcours du fichier)
        for (i,line) in enumerate(f):
            current_line = line.replace('\n','').split(delim)
            assert len(current_line) == len(new_type_list), "Error number of columns. Line:%d" % (i+1)
            yield tuple(t(x) for (x,t) in zip(current_line,new_type_list))
        f.close()

    @staticmethod
    def MySQLFileLoader(f):
        """load MySQL dumps (join truncated lines)"""
        tmp = ""
        for ligne in f:
            ligne = ligne.replace('\n', '')
            if ligne[-1] == '\\':
                # sign that shows that the line is not finished
                tmp = ligne[:-1]
            else:
                yield tmp + ligne
                tmp = ""
        assert (tmp == "")

    @staticmethod
    def MySQLFileWriter(data):
        """write a mySQL dump file"""
        # return myTSV.printLine(data).replace("None", "\N")
        # 2to3 (unicode error on \N)
        return myTSV.printLine(data).replace("None", "\\N")

class firstLineBuffer:
    """Read the first line of a file. Useful when you want to know the format
    without having to open-close it before opening it once more"""
    def __init__(self, f):
        self.f = f
        try:
            self.firstLine = next(self)
        except StopIteration:
            self.firstLine = ""

    def __iter__(self):
        yield self.firstLine
        while True:
            try:
                yield next(self)
            except StopIteration:
                return

    def __next__(self):
        while True:
            l = next(self.f).replace('\n', '').replace('\r', '')
            # Suppression of the lines with comments
            if (not l.startswith("#")) and (len(l) > 0):
                return l

    def next(self): # python2 compatibility
        return self.__next__()

    def close(self):
        return self.f.close()


def hasAccess(s):
    """existing file"""
    return os.access(os.path.expanduser(s), os.R_OK)


def openFile(nom, mode):
    """open a file and decompress it if possible
    return the object 'file' and the full name of the file"""

    # file already open
    if type(nom) != str:
        return nom

    # file on the web
    elif nom.startswith("http://") or nom.startswith("ftp://"):
        comm = "wget %s -O -"
        # Compression bzip2
        if nom.endswith(".bz2"):
            comm += " | bunzip2"
        # Compression gzip
        elif nom.endswith(".gz"):
            comm += " | gunzip"
        # Compression lzma
        elif nom.endswith(".lzma"):
            comm += " | unlzma"
        (stdin,f,stderr) = os.popen3( comm % nom )
        stdin.close()
        stderr.close()

    # standard entry
    elif nom == "-":
        return sys.stdin

    # file on the disk
    else:
        nom = os.path.expanduser(nom)
        if ("w" in mode) or ("a" in mode):
            # create the folder for the output into files
            try:
                os.makedirs(os.path.dirname(nom))
            except OSError:
                pass
        need_decoding = True
        i = nom.find(".zip/")
        if (mode == "r") and (i >= 0):
            import zipfile
            import io
            f = zipfile.ZipFile(nom[:i+4], "r")
            f = io.StringIO(f.read(nom[i+5:]))
        # Compression bzip2
        elif nom.endswith(".bz2"):
            import bz2
            f = bz2.BZ2File(nom, mode)
        # Compression gzip
        elif nom.endswith(".gz"):
            import gzip
            f = gzip.GzipFile(nom, mode)
        # Compression lzma
        elif nom.endswith(".lzma") or nom.endswith(".xz"):
            import lzma
            f = lzma.LZMAFile(nom, mode)
        else:
            f = open(nom, mode)
            need_decoding = False

        if need_decoding:
            f = codec_hook(f)
    return f
