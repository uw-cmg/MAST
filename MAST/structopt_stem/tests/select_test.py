from atompy.selection import *

def select_test():
	'''Selection function unit tests
	Key to test is:
		- results in anticipated order
		- no repetitive energies
	'''
	print 'Beginning unit testing of selection functions'
	try:
		print 'Best selection test SUCCESSFUL'
	except Exception, e:
		print 'ERROR: Best selection test FAILED'
		print e
		pass
	try:
		print 'Cost selection test SUCCESSFUL'
	except Exception, e::
		print 'ERROR: Cost selection test FAILED'
		print e
		pass
	try:
		print 'FUSS selection test SUCCESSFUL'
	except Exception, e:
		print 'ERROR: FUSS selection test FAILED'
		print e
		pass
	try:
		print 'FUSS1 selection test SUCCESSFUL'
	except Exception, e:
		print 'ERROR: FUSS1 selection test FAILED'
		print e
		pass
	try:
		print 'FUSS1R selection test SUCCESSFUL'
	except Exception, e:
		print 'ERROR: FUSS1R selection test FAILED'
		print e
		pass
	try:
		print 'FUSS2 selection test SUCCESSFUL'
	except Exception, e:
		print 'ERROR: FUSS2 selection test FAILED'
		print e
		pass
	try:
		print 'FUSSF selection test SUCCESSFUL'
	except Exception, e:
		print 'ERROR: FUSSF selection test FAILED'
		print e
		pass
	try:
		print 'FUSSR selection test SUCCESSFUL'
	except Exception, e:
		print 'ERROR: FUSSR selection test FAILED'
		print e
		pass
	try:
		print 'Metropolis selection test SUCCESSFUL'
	except Exception, e:
		print 'ERROR: Metropolis selection test FAILED'
		print e
		pass
	try:
		print 'MultiTournament selection test SUCCESSFUL'
	except Exception, e:
		print 'ERROR: MultiTournament selection test FAILED'
		print e
		pass
	try:
		print 'Random selection test SUCCESSFUL'
	except Exception, e:
		print 'ERROR: Random selection test FAILED'
		print e
		pass
	try:
		print 'Rank selection test SUCCESSFUL'
	except Exception, e:
		print 'ERROR: Rank selection test FAILED'
		print e
		pass
	try:
		print 'Tournament selection test SUCCESSFUL'
	except Exception, e:
		print 'ERROR: Tournament selection test FAILED'
		print e
		pass
	try:
		print 'Tournament1 selection test SUCCESSFUL'
	except Exception, e:
		print 'ERROR: Tournament1 selection test FAILED'
		print e
		pass
	try:
		print 'Tournament2 selection test SUCCESSFUL'
	except Exception, e:
		print 'ERROR: Tournament2 selection test FAILED'
		print e
		pass
	print 'Selection Function Testing Complete'
	return
	