from DAG import DAGScheduler as ds
myds = ds.DAGScheduler()

task_list = myds.schedule("//home/tam/bin/git/MAST4pymatgen/MAST/test/nebtest/test-recipe.txt")
print task_list
