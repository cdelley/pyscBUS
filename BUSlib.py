#
#   Author: Cyrille L. Delley
#   
#   Functions to read and write uncompressed BUS files with python.
#    
#   The file format was introduced by Páll Melsted, Sina Booeshaghi, Lior Pachtr
#   and coworkers for kallisto.
#
#   Melsted, Páll, Booeshaghi, A. Sina et al. Modular and efficient pre-processing
#   of single-cell RNA-seq. BioRxiv (2019): 673285, doi.org/10.1101/673285.   
#   See file specifications at:
#
#   https://github.com/BUStools/BUS-format
#   https://github.com/BUStools/bustools
#

import csv
import math
import struct
from typing import List, Dict, Tuple, Union, BinaryIO
import warnings

BUSFORMAT_VERSION = 1

def stringToBinary(s: str, length: int = 32, formating: str = '064b')  -> int:
    """
    Converts DNA string to binary number representation with:
    A   -> 00
    C   -> 01
    G   -> 10
    T   -> 11
    
    Input
        s: DNA string
        length: length of s, if not provided uses len(s) (default None)
        formating: binary format of the output, (default '032b')
    
    Returns
        int (formatted accoding to 'formating')
    """
    bs = s[:length].replace(
        'A','00').replace(
        'C','01').replace(
        'G','10').replace(
        'T','11')
    r = int(bs, 2)
    #return format(r, formating)
    return int(r)
    
def binaryToString(x: int, length: int)  -> str:
    """
    Reverts an binary integer number into a DNA string with:
    A   <- 00
    C   <- 01
    G   <- 10
    T   <- 11
    
    string length needs to be specified to get the appropriate number
    of leading 'A's
    """
    s = ["N"] * length
    sh = length - 1
    for i in range(length):
        c = "N"
        if ((x >> (2 * sh)) & 0x03) == 0x00:
            c = "A"
        elif ((x >> (2 * sh)) & 0x03) == 0x01:
            c = "C"
        elif ((x >> (2 * sh)) & 0x03) == 0x02:
            c = "G"
        elif ((x >> (2 * sh)) & 0x03) == 0x03:
            c = "T"
        sh -= 1
        s[i] = c
    return "".join(s)

class BUSFIle(object):
    def __init__(self):
        self.text = ""
        self.bcd = dict()  # type: dict[BC : [(UMI1, ecs1, count, flag), 
                          #                  (UMI2, ecs1, count, flag), 
                          #                  (UMI1, ecs2, count, flag), ....]
        self.version = 0
        self.bclen = 0
        self.umilen = 0
        self.textlen = 0
        self.n_entries = 0
        self.format = ''
    
    @classmethod
    def read(cls, in_file: str) -> "BUSFIle":
        """
        read BUS uncompressed files and return a BUSFile instance
        """
        instance = cls()
        instance.in_file_path = in_file
        with open(in_file, 'rb') as inf:
            instance.parseHeader(inf)
            if instance.format == 'BUS':
                instance.parseBUS(inf)
            else:
                print('file format {} not implemented'.format(instance.format))
        return instance

    def parseBUS(self, inf: BinaryIO) -> None:
        """
        parse the data blocks of an uncompressed BUS file
        """
        b = inf.read(32)
        self.n_entries += 1
        while b:
            bc, umi, ec, num, flag, pad = struct.unpack('<QQIIII', b)
            b = inf.read(32)
            self.n_entries += 1
            try:
                self.bcd[bc].append((umi, ec, num, flag))
            except KeyError:
                self.bcd[bc] = [(umi, ec, num, flag)] 
    
    def BUS_to_txt(self) -> None:
        # binaryToString(bc, self.bclen), binaryToString(umi, self.umilen), ec, num
        pass
        
    def parseHeader(self, inf: BinaryIO) -> None:
        """
        parse the header of an uncompressed BUS file
        """
        magic = struct.unpack("<4s", inf.read(4))[0]
        if magic == b"BUS\x00":
            self.format = 'BUS'
        else:
            self.format = magic    
            return       
        self.version, self.bclen, self.umilen, self.textlen = struct.unpack('<IIII', inf.read(16))
        if self.version != BUSFORMAT_VERSION:
            warnings.warn("The BUS file version is different from the PyBUS decoder version")
        self.text = struct.unpack("<{}s".format(self.textlen), inf.read(self.textlen))[0]
    
    def write(self, outf: str, formating: str = 'BUS') -> None:
        """
        write an uncompressed BUS file
        """
        self.out_file_path = outf
        with open(outf, 'wb') as fout:
            self.writeHeader(fout)
            if formating == 'BUS':
                self.writeBUS(fout)
    
    def writeHeader(self, outf: BinaryIO) -> None:
        """
        write the header of an uncompressed BUS file
        """
        try:
            text_bytes = self.text.encode()
        except AttributeError:
            text_bytes = self.text
        outf.write(b"BUS\x00")
        outf.write(struct.pack('<IIII', self.version, self.bclen, self.umilen, len(text_bytes)))
        outf.write(text_bytes)
        
    def writeBUS(self, outf: BinaryIO, sort_key: str = None) -> None:
        """
        write the data blocks of an uncompressed BUS file
        """
        if not sort_key:
            for bc, value in self.bcd.items():
                for tup in value:
                    block = struct.pack('<QQIIII', bc, tup[0], tup[1], tup[2], tup[3], 0)
                    outf.write(block)
            

