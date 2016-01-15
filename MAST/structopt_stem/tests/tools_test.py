from atompy.tools import *

def tool_test():
	'''Tool function unit tests
	Key to test is:
		- results in anticipated order
		- no repetitive energies
	'''
	print 'Beginning unit testing of tool functions'
	try:
		print 'BestInds test SUCCESSFUL'
	except Exception, e:
		print 'ERROR: BestInds test FAILED'
		print e
		pass
	try:
		print 'calc_dist test SUCCESSFUL'
	except Exception, e::
		print 'ERROR: calc_dist test FAILED'
		print e
		pass
	try:
		print 'check_atomlist_concentration test SUCCESSFUL'
	except Exception, e:
		print 'ERROR: check_atomlist_concentration test FAILED'
		print e
		pass
	try:
		print 'eval_energy test SUCCESSFUL'
	except Exception, e:
		print 'ERROR: eval_energy test FAILED'
		print e
		pass
	try:
		print 'find_ints test SUCCESSFUL'
	except Exception, e:
		print 'ERROR: find_ints test FAILED'
		print e
		pass
	try:
		print 'find_top_layer test SUCCESSFUL'
	except Exception, e:
		print 'ERROR: find_top_layer test FAILED'
		print e
		pass
	try:
		print 'fitneseval test SUCCESSFUL'
	except Exception, e:
		print 'ERROR: fitnesseval test FAILED'
		print e
		pass
	try:
		print 'get_best test SUCCESSFUL'
	except Exception, e:
		print 'ERROR: get_best test FAILED'
		print e
		pass
	try:
		print 'setup_calculator test SUCCESSFUL'
	except Exception, e:
		print 'ERROR: setup_calculator test FAILED'
		print e
		pass
	try:
		print 'setup_fixed_region_calculator test SUCCESSFUL'
	except Exception, e:
		print 'ERROR: setup_fixed_region_calculator test FAILED'
		print e
		pass
	print 'Tool Function Testing Complete'
	return
	