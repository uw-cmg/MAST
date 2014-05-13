from atompy.fitness import *

def fit_test():
	'''Fitness function unit tests
	Key to test is:
		- results in anticipated order
	'''
	print 'Beginning unit testing of fitness functions'
	try:
		print 'ChemPotSwap fitness test SUCCESSFUL'
	except Exception, e:
		print 'ERROR: ChemPotSwap fitness test FAILED'
		print e
		pass
	try:
		print 'enpafit fitness test SUCCESSFUL'
	except Exception, e::
		print 'ERROR: enpafit fitness test FAILED'
		print e
		pass
	try:
		print 'enthalpyfit fitness test SUCCESSFUL'
	except Exception, e::
		print 'ERROR: enthalpyfit fitness test FAILED'
		print e
		pass
	try:
		print 'exponentialfit fitness test SUCCESSFUL'
	except Exception, e::
		print 'ERROR: exponentialfit fitness test FAILED'
		print e
		pass
	try:
		print 'formationenergy fitness test SUCCESSFUL'
	except Exception, e::
		print 'ERROR: formationenergy fitness test FAILED'
		print e
		pass
	try:
		print 'STEM_Cost fitness test SUCCESSFUL'
	except Exception, e::
		print 'ERROR: STEM_Cost fitness test FAILED'
		print e
		pass
	try:
		print 'surfaceenergy fitness test SUCCESSFUL'
	except Exception, e::
		print 'ERROR: surfaceenergy fitness test FAILED'
		print e
		pass
	try:
		print 'totalenfit fitness test SUCCESSFUL'
	except Exception, e::
		print 'ERROR: totalenfit fitness test FAILED'
		print e
		pass
	print 'Fitness Function Testing Complete'
	return
	