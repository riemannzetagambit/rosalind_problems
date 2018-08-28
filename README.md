**Overview**
All problems from Rosalind are at a given URL `http://rosalind.info/problems/{PROBLEM_NAME}/`.
The solution to each problem can be found at the corresponding python file `rosalind_{PROBLEM_NAME}.py`.
For clarity, I do not include the associated `.txt` files of the same basename and only store them locally.

In addition, there is a single module `rosalind_utils` that I include that has convenience functions and data that I use
widely across problems-- e.g. the D/RNA alphabets, translation dictionaries, and some helper functions for testing.

**Getting problem solutions**
All of the solution scripts `rosalind_{PROBLEM_NAME}.py` are set up to take a single command line argument-- the data
file provided by Rosalind.

A sample usage for the very first problem would be
```python
python rosalind_dna.py rosalind_dna.txt
```
All outputs should give you solutions in the format that Rosalind requests.

**Compatibility**
I assume python 3 and the standard library for all solutions.
