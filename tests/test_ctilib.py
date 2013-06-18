from unittest import TestCase, main
from os import tmpfile, remove
from multiprocessing import Manager, Pool
from selextrace.ctilib import fold_clusters, group_by_shape, \
    run_fold_for_infernal, build_reference, group_to_reference, group_denovo, \
    group_by_forester, make_r2r, run_infernal


class MainTests(TestCase):
    def setUp(self):
        self.fastaseqs = {"MISEQ:8:000000000-A18AY:1:1101:12274:26426_11":
            "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCGUGGAUAUGGUACACAUUGAAUUGUUGGACACCGUAAAUGUCCUAACACGGGGGCAU",
            "MISEQ:8:000000000-A18AY:1:1101:14188:4443_26":
            "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCGUGGAUAUGGCGCGCAUUGAAUUGUUGGGCACCGUAAAUGUCCUAACACGGGGGCAU",
            "MISEQ:8:000000000-A18AY:1:1101:22090:7020_90":
            "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCGUGGAUAUGGUACGCAUUGAAUUGUUGGGCACCGUAAAUGUCCUAACACGGGGGCAU",
            "MISEQ:8:000000000-A18AY:1:1102:15083:7079_17":
            "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCAUGGAUAUGGUACGCAUUGAAUUGUUGGGCACCGUAAAUGUCCUAACACGGGCGCAU",
            "MISEQ:8:000000000-A18AY:1:1102:18612:5535_7":
            "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCGUGGAUAUGGUACGCAUUGAAUUGUUGGGCACCGUAAAUGCCCUAACACGGGGGCAU",
            "MISEQ:8:000000000-A18AY:1:1102:26024:20405_4":
            "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCGUGGAUAUGGCCCGCAUUGAAUUGUUGGGCACCGUAAAUGUCCUAACACGGGGGCAU",
            "MISEQ:8:000000000-A18AY:1:1105:14627:21327_3":
            "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCAUGGAUAUGGUACGCAUUGAAUUGUUGGGCACCGUAAAUGCCCUAACACGGGCGCAU",
            "MISEQ:8:000000000-A18AY:1:1105:9939:17358_2":
            "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCAUGGAUAUGGUACGCAUUGAAUUGUUGGGCACCGUAAAUGUCCUAACACGGGCACAU",
            "MISEQ:8:000000000-A18AY:1:1106:17400:6900_2":
            "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGAGUGGAUAUGGUACGCAUUGAAUUGUUGGGCACCGUAAAUGUCCUAACACGGGGGCAU",
            "MISEQ:8:000000000-A18AY:1:1107:2414:17210_2":
            "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAGCGCGUGGAUAUGGUACGCAUUGAAUUGUUGGACACCGUAAAUGUCCUAACACGGGGGCAU",
            "MISEQ:8:000000000-A18AY:1:1107:29080:16480_2":
            "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCGUGGAUAUGGUACGCAUUGAAUUGUUGGGCACCGUAAAUGUCCUAACACGGGGGCAA",
            "MISEQ:8:000000000-A18AY:1:1110:11594:5021_2":
            "GACUUCGGUCCAAGCUAAUGCACUCUGAUGAUCGCGUGGAUAUGGUACGCAUUGAAUUGUUGGGCACCGUAAAUGUCCUAACCCGGGGGCAU"}
        self.fastastruct = ".....((((((..((((((((...((....(((.....)))..))...))))))....))..)).))))....((((((......))))))."

    def test_fold_clusters(self):
        cfo = open("cluster_structs.fasta", 'w')
        cfo.close()
        manager = Manager()
        lock = manager.Lock()
        fold_clusters(lock, "cluster_1", self.fastaseqs, "./")
        obs = ''.join(open("cluster_structs.fasta").read())
        remove("cluster_structs.fasta")
        expected = ">cluster_1\n" + self.fastastruct + "\n"
        self.assertEqual(obs, expected)

    def test_group_by_shape_add(self):
        #run the pool over all groups to get structures
        manager = Manager()
        obs = manager.dict()
        obs["[]"] = ["(((...)))"]
        pool = Pool(processes=1)
        pool.apply_async(func=group_by_shape,
            args=(obs, "(((((((((((..........)))))))))))"))
        pool.close()
        pool.join()
        exp = {
            "[]": ["(((...)))", "(((((((((((..........)))))))))))"]
        }
        for shape in exp:
            if shape not in obs.keys():
                raise AssertionError(shape + " not in observed!")
            self.assertEqual(obs[shape], exp[shape])

    def test_group_by_shape_new(self):
        #shape [][]
        manager = Manager()
        obs = manager.dict()
        obs["[]"] = ["(((...)))"]
        pool = Pool(processes=1)
        pool.apply_async(func=group_by_shape,
            args=(obs, "..(((.....))).....(((...))).."))
        pool.close()
        pool.join() 
        exp = {
            "[]": ["(((...)))"],
            "[][]": ["..(((.....))).....(((...))).."]
        }
        for shape in exp:
            if shape not in obs.keys():
                raise AssertionError(shape + " not in observed!")
            self.assertEqual(obs[shape], exp[shape])


    #def test_build_reference(self):


    #def test_group_to_reference(self):


    #def test_group_denovo(self):


    #def test_run_for_infernal(self):


    #def test_group_by_forester(self):


    #def test_make_r2r(self):


    #def test_run_infernal(self):


if __name__ == "__main__":
    main()