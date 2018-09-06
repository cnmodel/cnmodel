# -*- encoding: utf-8 -*-

import sys
import os
import unittest


class Params(object):
    def __init__(self, **kwds):
        """
        Utility class to create parameter lists
        create using:
        p = Params(abc=2.0, defg = 3.0, lunch='sandwich')
        reference using:
        p.abc, p.defg, etc.
        Supports getting the keys, finding whether a key exists, returning the strucure as a simple dictionary,
        and printing (show) the parameter structure.
        """
        self.__dict__.update(kwds)

    def additem(self, key, value):
        self.__dict__[key] = value

    def getkeys(self):
        """
        Get the keys in the current dictionary
        """
        return(self.__dict__.keys())

    def haskey(self, key):
        """
        Find out if the param list has a specific key in it
        """
        if key in self.__dict__.keys():
            return True
        else:
            return False

    def todict(self):
        """
        convert param list to standard dictionary
        Useful when writing the data
        """
        r = {}
        for dictelement in self.__dict__:
            if isinstance(self.__dict__[dictelement], Params):
                #print 'nested: ', dictelement
                r[dictelement] = self.__dict__[dictelement].todict()
            else:
                r[dictelement] = self.__dict__[dictelement]
        return r

    def show(self, printFlag = True):
        """
        print the parameter block created in Parameter Init
        """
        print("--------    Parameter Block    ----------")
        for key in self.__dict__.keys():
            print("%15s = " % (key), eval('self.%s' % key))
        print("-------- ---------------------- ----------")


class ParamTests(unittest.TestCase):
	def setUp(self):
		pass


if __name__ == '__main__':
	unittest.main()