#!/usr/bin/env python
class TestClass():
    def __init__(self, name="banana"):
        print "Hello!"
        self.name=name
    def run_me(self):
        print "Test ran with value %s" % self.name
import sys
print "System argument: ", sys.argv
myclass = TestClass("pineapple")
callmethod = getattr(TestClass, sys.argv[1])
result=callmethod(myclass)
result=getattr(TestClass, sys.argv[1])(myclass)
#import foo
#methodToCall = getattr(foo, 'bar')
#result = methodToCall()
#
#As far as that goes, lines 2 and 3 can be compressed to:
#
#result = getattr(foo, 'bar')()

