from unittest import TestCase, main
from selextrace.remove_duplicates import remove_duplicates


class MainTests(TestCase):
        def test_remove(self):
            testin = [
                ("test1", "AAAAGGGCCCTTTAGCTAAA"),
                ("test2", "AAAGGGCCCTTTAAA"),
                ("test3", "AAAGGGCCCTTTAAA")]
            obs1, obs2 = remove_duplicates(testin)

            expected = [('test2_2', 'AAAGGGCCCTTTAAA'),
                ('test1_1', 'AAAAGGGCCCTTTAGCTAAA')]
            self.assertEqual(obs1, expected)

        def test_remove2(self):
            testin = [
                ("test1", "AAAAGGGCCCTTTAGCTAAA"),
                ("test2", "AAAGGGCCCTTTAAA"),
                ("test3", "AAAGGGCCCTTTAAA")]
            obs1, obs2 = remove_duplicates(testin)
            expected = {
                "AAAAGGGCCCTTTAGCTAAA": ["test1"],
                "AAAGGGCCCTTTAAA": ["test2", "test3"]}
            self.assertEqual(obs2, expected)


if __name__ == "__main__":
    main()
