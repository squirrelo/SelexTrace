from unittest import TestCase, main
from remove_duplicates import remove_duplicates

class MainTests(TestCase):
        def test_remove(self):
            testin = [("test1","AAAAGGGCCCTTTAGCTAAA"), ("test2","AAAGGGCCCTTTAAA"), ("test3","AAAGGGCCCTTTAAA")]
            test1, test2= remove_duplicates(testin)
            self.assertEqual(test1,[('test2_2', 'AAAGGGCCCTTTAAA'), ('test1_1', 'AAAAGGGCCCTTTAGCTAAA')])
        def test_remove2(self):
            testin = [("test1","AAAAGGGCCCTTTAGCTAAA"), ("test2","AAAGGGCCCTTTAAA"), ("test3","AAAGGGCCCTTTAAA")]
            test1, test2= remove_duplicates(testin)
            self.assertEqual(test2,{"AAAAGGGCCCTTTAGCTAAA":["test1"], "AAAGGGCCCTTTAAA":["test2", "test3"]})


if __name__ == "__main__":
    main()
