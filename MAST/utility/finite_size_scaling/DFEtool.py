from MAST.utility.finite_size_scaling.defectformationenergy import DefectFormationEnergy
from MAST.utility.finite_size_scaling.cubicscaling import CubicScaling
from MAST.utility.finite_size_scaling.defectlevels import DefectLevels
from MAST.utility.finite_size_scaling.parse_input import ParseInput

if __name__=="__main__":
    inp = ParseInput().read_input()
    GGA = inp['GGA']
    HSE = inp['HSE']
    method = inp['method']
    DefectFormationEnergy().main(GGA,HSE)
    for chempot in GGA['chem_pot'].keys(): 
        CubicScaling(inputfile='%s_%s.txt'%(chempot,GGA['tag']),method=method).main()
        DefectLevels(inputfile='%s_%s.txt'%(chempot,method),XC=[GGA,HSE]).main()
    import os,shutil
    os.mkdir('dfe_results')
    for file in os.listdir('.'):
        if (not 'metadata' in file) and (not 'inp' in file) and ('.txt' in file) or ('.eps' in file) or ('.stat' in file):
            shutil.move(file,'dfe_results')

