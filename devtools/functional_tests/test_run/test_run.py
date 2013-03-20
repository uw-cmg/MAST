"""Tests for Run"""

from MAST.control.run import Run
from MAST.interface.wholeinput import WholeInput
from MAST.mastfile.mastfile import MASTFile
from MAST.mastfile.kpoints import Kpoints
from MAST.mastfile.poscar import Poscar
from MAST.mastfile.potcar import Potcar
from MAST.mastfile.incar import Incar
import os
import time
import subprocess
from nose.plugins.skip import SkipTest

class TestRun:

    def setUp(self):
        stem=os.getcwd()
        self.finished_path = stem + "/finished_run"
        self.unfinished_path = stem + "/unfinished_run"
        self.unfinishednc_path = stem + "/unfinished_run_nocontcar"
        self.unfinishedec_path = stem + "/unfinished_run_emptycontcar"
        self.nopotcar_path = stem + "/unfinished_run_nopotcar"
        self.noposcar_path = stem + "/unfinished_run_noposcar"
        self.nokpoints_path = stem + "/unfinished_run_nokpoints"
        self.noincar_path = stem + "/unfinished_run_noincar"
        self.noscript_path = stem + "/unfinished_run_noscript"
        self.test_cgin = WholeInput()

    def tearDown(self):
        try:
            os.remove(self.unfinished_path + "/POSCAR1")
            os.remove(self.unfinished_path + "/CONTCAR1")
            os.remove(self.unfinished_path + "/OUTCAR1")
            os.remove(self.unfinishednc_path + "/OUTCAR1")
        except OSError:
            pass
    
    def test___init__(self):
        self.test_cgin["calcgroup_path"] = self.finished_path
        self.testclass = Run(self.test_cgin)
        assert self.testclass.status == 'I'
        assert self.testclass.cgin == self.test_cgin

    def test_create_essential_files(self):
        self.test_cgin["calcgroup_path"] = self.finished_path
        self.testclass = Run(self.test_cgin)
        self.testclass.runpath = self.test_cgin["calcgroup_path"]
        #
        self.testclass.create_essential_files()
        #
        testk = Kpoints(self.testclass.runpath + "/KPOINTS")
        assert testk.data == self.testclass.kpoints.data
        testi = Incar(self.testclass.runpath + "/INCAR")
        testi.create_data() #necessary for INCAR
        assert testi.data == self.testclass.incar.data
        testpt = Potcar(self.testclass.runpath + "/POTCAR")
        assert testpt.data == self.testclass.potcar.data
        testpc = Poscar(self.testclass.runpath + "/POSCAR")
        assert testpc.data == self.testclass.poscar.data

    def test_set_up_output_paths(self):
        self.test_cgin["calcgroup_path"] = self.finished_path
        self.testclass=Run(self.test_cgin)
        self.testclass.runpath = self.test_cgin["calcgroup_path"]
        #
        self.testclass.set_up_output_paths()
        #   
        assert self.testclass.contcar.path == self.testclass.runpath + "/CONTCAR"
        assert self.testclass.outcar.path == self.testclass.runpath + "/OUTCAR"
        assert self.testclass.oszicar.path == self.testclass.runpath + "/OSZICAR"

    def test_clean_directory(self):
        raise SkipTest
        #self.testclass.clean_directory()

    def test_missing_restart_files(self):
        #No files missing
        self.test_cgin["calcgroup_path"] = self.finished_path
        self.testclass=Run(self.test_cgin)
        self.testclass.runpath = self.test_cgin["calcgroup_path"]
        self.testclass.set_up_output_paths()
        assert self.testclass.missing_restart_files() == None
        #CONTCAR missing
        self.test_cgin["calcgroup_path"] = self.unfinishednc_path
        self.testclass=Run(self.test_cgin)
        self.testclass.runpath = self.test_cgin["calcgroup_path"]
        self.testclass.set_up_output_paths()
        assert self.testclass.missing_restart_files() == True
        #CONTCAR empty
        self.test_cgin["calcgroup_path"] = self.unfinishedec_path
        self.testclass=Run(self.test_cgin)
        self.testclass.runpath = self.test_cgin["calcgroup_path"]
        self.testclass.set_up_output_paths()
        assert self.testclass.missing_restart_files() == True


    def test_missing_essential_files(self):
        #No files missing
        self.test_cgin["calcgroup_path"] = self.finished_path
        self.testclass=Run(self.test_cgin)
        self.testclass.runpath = self.test_cgin["calcgroup_path"]
        assert self.testclass.missing_essential_files() == None
        #Missing INCAR
        self.test_cgin["calcgroup_path"] = self.noincar_path
        self.testclass=Run(self.test_cgin)
        self.testclass.runpath = self.test_cgin["calcgroup_path"]
        assert self.testclass.missing_essential_files() == True
        #Missing POTCAR
        self.test_cgin["calcgroup_path"] = self.nopotcar_path
        self.testclass=Run(self.test_cgin)
        self.testclass.runpath = self.test_cgin["calcgroup_path"]
        assert self.testclass.missing_essential_files() == True
        #Missing POSCAR
        self.test_cgin["calcgroup_path"] = self.noposcar_path
        self.testclass=Run(self.test_cgin)
        self.testclass.runpath = self.test_cgin["calcgroup_path"]
        assert self.testclass.missing_essential_files() == True
        #Missing KPOINTS
        self.test_cgin["calcgroup_path"] = self.nokpoints_path
        self.testclass=Run(self.test_cgin)
        self.testclass.runpath = self.test_cgin["calcgroup_path"]
        assert self.testclass.missing_essential_files() == True


    def test_missing_submission_script(self):
        raise SkipTest
        #self.testclass.missing_submission_script()

    def test_find_last_submission_script(self):
        self.test_cgin["calcgroup_path"] = self.finished_path
        self.testclass=Run(self.test_cgin)
        self.testclass.runpath = self.test_cgin["calcgroup_path"]
        assert self.testclass.find_last_submission_script() != None
        #no script
        self.test_cgin["calcgroup_path"] = self.noscript_path
        self.testclass=Run(self.test_cgin)
        self.testclass.runpath = self.test_cgin["calcgroup_path"]
        assert self.testclass.find_last_submission_script() == None


    def test_okay_to_submit(self):
        #
        self.test_cgin["calcgroup_path"] = self.finished_path
        self.testclass=Run(self.test_cgin)
        self.testclass.runpath = self.test_cgin["calcgroup_path"]
        assert self.testclass.okay_to_submit() == False
        #
        self.test_cgin["calcgroup_path"] = self.unfinished_path
        self.testclass=Run(self.test_cgin)
        self.testclass.runpath = self.test_cgin["calcgroup_path"]
        assert self.testclass.okay_to_submit() == True
        #
        self.test_cgin["calcgroup_path"] = self.unfinishednc_path
        self.testclass=Run(self.test_cgin)
        self.testclass.runpath = self.test_cgin["calcgroup_path"]
        assert self.testclass.okay_to_submit() == False
        #

    def test_queue_status(self):
        raise SkipTest
        #self.testclass.queue_status(self.testclass.cgin, jobid)

    def test_submit_to_queue(self):
        raise SkipTest
        #self.testclass.submit_to_queue()

    def test_store_jobid(self):
        raise SkipTest
        #self.testclass.store_jobid()

    def test_get_latest_jobid(self):
        raise SkipTest
        #self.testclass.get_latest_jobid()

    def test_recycle_contcar(self):
        raise SkipTest
        #self.testclass.recycle_contcar()

    def test_archive_output(self):
        raise SkipTest
        #self.testclass.archive_output()

    def test_updated_within_lagseconds(self):
        raise SkipTest
        #self.testclass.updated_within_lagseconds()

    def test_record_to_runfile(self):
        raise SkipTest
        #self.testclass.record_to_runfile()

    def test_update_status(self):
        raise SkipTest
        #self.testclass.update_status()

    def test_parse_efile(self):
        raise SkipTest
        #self.testclass.parse_efile()

    def test_is_complete(self):
        #
        self.test_cgin["calcgroup_path"] = self.finished_path
        self.testclass=Run(self.test_cgin)
        self.testclass.runpath = self.test_cgin["calcgroup_path"]
        assert self.testclass.is_complete() == True
        #
        self.test_cgin["calcgroup_path"] = self.unfinished_path
        self.testclass=Run(self.test_cgin)
        self.testclass.runpath = self.test_cgin["calcgroup_path"]
        assert self.testclass.is_complete() == False

    def test_print_finished(self):
        raise SkipTest
        #self.testclass.print_finished()

    def test_get_and_set_energy(self):
        self.test_cgin["calcgroup_path"] = self.finished_path
        self.testclass=Run(self.test_cgin)
        self.testclass.runpath = self.test_cgin["calcgroup_path"]
        self.testclass.set_up_output_paths()
        self.testclass.get_and_set_energy()
        assert self.testclass.energy == float('-.34198451E+03')

    def test_get_and_set_volume(self):
        self.test_cgin["calcgroup_path"] = self.finished_path
        self.testclass=Run(self.test_cgin)
        self.testclass.runpath = self.test_cgin["calcgroup_path"]
        self.testclass.set_up_output_paths()
        self.testclass.get_and_set_volume()
        assert self.testclass.volume == float(500.63)

    def test_get_and_set_num_atoms(self):
        self.test_cgin["calcgroup_path"] = self.finished_path
        self.testclass=Run(self.test_cgin)
        self.testclass.runpath = self.test_cgin["calcgroup_path"]
        self.testclass.create_essential_files()
        self.testclass.set_up_output_paths()
        self.testclass.get_and_set_num_atoms()
        assert self.testclass.num_atoms == 39

    def test_get_and_set_energy_per_atom(self):
        self.test_cgin["calcgroup_path"] = self.finished_path
        self.testclass=Run(self.test_cgin)
        self.testclass.runpath = self.test_cgin["calcgroup_path"]
        self.testclass.create_essential_files()
        self.testclass.set_up_output_paths()
        self.testclass.get_and_set_energy_per_atom()
        assert self.testclass.energy == float('-.34198451E+03')
        assert self.testclass.energy_per_atom == float('-.34198451E+03')/39

    def test_get_and_set_volume_per_atom(self):
        self.test_cgin["calcgroup_path"] = self.finished_path
        self.testclass=Run(self.test_cgin)
        self.testclass.runpath = self.test_cgin["calcgroup_path"]
        self.testclass.set_up_output_paths()
        self.testclass.create_essential_files()
        self.testclass.get_and_set_volume_per_atom()
        assert self.testclass.volume == float(500.63)
        assert self.testclass.volume_per_atom == float(500.63)/39

