from MAST.DAG.jobentry import JobEntry

def test_jobentry_init():
    je = JobEntry(1,2,jobname='/mast/scratch/recipe/eachjob/',type=None, ingredient=None)
    print je

def test_jobentry_methods():
    je = JobEntry(1,2,jobname='/mast/scratch/recipe/eachjob/',type=None, ingredient=None)
    print je
    je.addparent(3)
    je.addparent(4)
    je.addparent(5)
    print je 
    assert je.is_ready() == False
    je.completeparent(3)
    je.completeparent(4)
    je.completeparent(5)
    print je
    assert je.is_ready() == True
    
if __name__ == "__main__":
    test_jobentry_init()
    test_jobentry_methods()
    
