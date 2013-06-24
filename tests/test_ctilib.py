from unittest import TestCase, main
from os import tmpfile, remove
from multiprocessing import Manager, Pool
from selextrace.ctilib import fold_clusters, group_by_shape, \
    run_fold_for_infernal, build_reference, group_to_reference, group_denovo, \
    score_local_rnaforester, group_by_forester, group_by_distance, make_r2r, run_infernal


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
        self.testclusterslt10 = {
            "(((...)))": [(1, 2), (3, 4)],
            "((.....))": [(5, 6)],
            ".....((((.......))))...": [(7, 8)],
            "(((...)))...(((...)))": [(9, 10)],
            "((((....))))": [(11, 12)],
            '............': [(13, 14)],
            '....((((((.......))))))....': [(15, 16)]
        }

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

    def test_build_reference_lengths(self):
        items = [1, 2, 3, 4, 5, 6, 7, 8, 9, 0]
        obs1, obs2 = build_reference(items, 3)
        self.assertEqual(len(obs1), 3)
        self.assertEqual(len(obs2), 7)

    def test_build_reference_items(self):
        items = [1, 2, 3, 4, 5, 6, 7, 8, 9, 0]
        obs1, obs2 = build_reference(items, 3)
        self.assertEqual(len(obs1), 3)
        self.assertEqual(len(obs2), 7)
        #check that the ref and nonref are same
        for num in obs1:
            if num not in items:
                raise AssertionError(str(num) + " not in expected!")
        for num in obs2:
            if num not in items:
                raise AssertionError(str(num) + " not in expected!")
        for num in items:
            if num not in obs1 and num not in obs2:
                raise AssertionError(str(num) + " not in observed!")

    def test_build_reference_dupes(self):
        items = [1, 2, 3, 4, 5, 6, 7, 8, 9, 0]
        obs1, obs2 = build_reference(items, 3)
        self.assertEqual(len(obs1), 3)
        self.assertEqual(len(obs2), 7)
        #check that the ref and nonref are same
        finals = set([])
        for item in obs1:
            if item in finals:
                raise AssertionError("Duplicate in reference!")
            finals.add(item)
        for item in obs2:
            if item in finals:
                raise AssertionError("Duplicate in nonreference!")
            finals.add(item)        

    def test_group_to_reference(self):
        testref = ["(((...)))"]
        testnonref = ["((.....))",
            ".....((((.......))))...",
            "(((...)))...(((...)))",
            "((((....))))",
            '............',
            '....((((((.......))))))....']
        exp1 = {'(((...)))...(((...)))': [(9, 10)],
        '(((...)))': [(1, 2), (3, 4), (5, 6), (11, 12)],
        '.....((((.......))))...': [(7, 8)],
        '............': [(13, 14)]}
        exp2 = ['.....((((.......))))...', '(((...)))...(((...)))', '............']
        obs1, obs2 = group_to_reference(self.testclusterslt10, testref, testnonref, 10)
        self.assertEqual(obs1, exp1)
        self.assertEqual(obs2, exp2)

    def test_group_denovo_all(self):
        obs1, obs2 = group_denovo(self.testclusterslt10, self.testclusterslt10.keys(), 10)
        print "\n", obs1
        print obs2
        exp1 = {
            '(((...)))...(((...)))': [(9, 10)],
            '(((...)))': [(1, 2), (3, 4), (11, 12), (5, 6)],
            '.....((((.......))))...': [(7, 8)],
            '............': [(13, 14)]}
        exp2 = ['.....((((.......))))...', '(((...)))...(((...)))', '(((...)))', '............']
        self.assertEqual(obs1, exp1)
        self.assertEqual(obs2, exp2)

    def test_group_denovo_specific(self):
        obs1, obs2 = group_denovo(self.testclusterslt10, ['.....((((.......))))...'], 10)
        exp1 = {
            '(((...)))...(((...)))': [(9, 10)],
            '(((...)))': [(1, 2), (3, 4), (11, 12), (5, 6)],
            '.....((((.......))))...': [(15, 16), (7, 8)],
            '............': [(13, 14)]}
        exp2 = ['.....((((.......))))...', '(((...)))...(((...)))', '(((...)))', '............']
        self.assertEqual(obs1, exp1)
        self.assertEqual(obs2, exp2)


    def test_group_by_distance_all(self):
        obs = group_by_distance(self.testclusterslt10, 10)
        exp = {
            '(((...)))...(((...)))': [(9, 10)],
            '(((...)))': [(1, 2), (3, 4), (11, 12), (5, 6)],
            '.....((((.......))))...': [(15, 16), (7, 8)],
            '............': [(13, 14)]}
        self.assertEqual(obs, exp)

    def test_group_by_distance_lt10_specific_nogroup(self):
        obs = group_by_distance(self.testclusterslt10, 10, ['.....((((.......))))...'])
        exp = {
            "(((...)))": [(1, 2), (3, 4)],
            "((.....))": [(5, 6)],
            ".....((((.......))))...": [(7, 8)],
            "(((...)))...(((...)))": [(9, 10)],
            "((((....))))": [(11, 12)],
            '............': [(13, 14)]
        }
        self.assertEqual(obs, exp)

    def test_group_by_distance_lt10_specific_group(self):
        obs = group_by_distance(self.testclusterslt10, 15, specstructs=['(((...)))'])
        print obs

    def test_score_local_rnaforester(self):
        struct1 = "....((.((.(((((((.((((((((....((...))(((((((...(((...((......))...))).))))))).....)))))).........))..))))).)).))))..."
        struct2 = "(((....))).......(((((..(((((......(((((((...(((...(((....)))...))).)))))))...)))))))))).......(((((((....))))))).."
        obs = score_local_rnaforester(struct1, struct2)
        self.assertEqual(obs, 137)


    #def test_group_by_rnaforester(self):

    #def test_run_for_infernal(self):

    #def test_make_r2r(self):

    #def test_run_infernal(self):


if __name__ == "__main__":
    main()