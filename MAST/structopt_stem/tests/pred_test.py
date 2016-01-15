from atompy.predator import *

def predator_test():
	'''Predator function unit tests
	Key to test is:
		- results in anticipated order
		- no repetitive energies
	'''
	print 'Beginning unit testing of predator functions'
	try:
		print 'Adapting predator test SUCCESSFUL'
	except Exception, e:
		print 'ERROR: Adapting predator test FAILED'
		print e
		pass
	try:
		print 'EnergyCluster predator test SUCCESSFUL'
	except Exception, e::
		print 'ERROR: EnergyCluster predator test FAILED'
		print e
		pass
	try:
		print 'Fingerprint_Niche predator test SUCCESSFUL'
	except Exception, e:
		print 'ERROR: Fingerprint_Niche predator test FAILED'
		print e
		pass
	try:
		print 'Fitpred_mut predator test SUCCESSFUL'
	except Exception, e:
		print 'ERROR: Fitpred_mut predator test FAILED'
		print e
		pass
	try:
		print 'Fitpred_new predator test SUCCESSFUL'
	except Exception, e:
		print 'ERROR: Fitpred_new predator test FAILED'
		print e
		pass
	try:
		print 'Fitpred predator test SUCCESSFUL'
	except Exception, e:
		print 'ERROR: Fitpred predator test FAILED'
		print e
		pass
	try:
		print 'Mutation_Dups_Energy predator test SUCCESSFUL'
	except Exception, e:
		print 'ERROR: Mutation_Dups_Energy predator test FAILED'
		print e
		pass
	try:
		print 'Mutation_Dups predator test SUCCESSFUL'
	except Exception, e::
		print 'ERROR: Mutation_Dups predator test FAILED'
		print e
		pass
	try:
		print 'ZPMutation_Dups predator test SUCCESSFUL'
	except Exception, e::
		print 'ERROR: ZPMutation_Dups predator test FAILED'
		print e
		pass
	print 'predator Function Testing Complete'
	return
	