# Wordpacker
Given a shape like a wordsquare, finds words that fit that shape.

You provide it with a topology file and a dictionary (one word per line, [a-z] only), Some example topology files:

Search for wordsquares size 4: (eg: each row represents a word, each number a letter, same numbers must be same letters)
```
1,2,3,4
2,5,6,7
3,6,8,9
4,7,9,10
```

Double word square of size 3
```
1,2,3
4,5,6
7,8,9
1,4,7
2,5,8
3,6,9
```

SATOR-like palinromic square:
```
1,2,3,4,5
2,6,7,8,4
3,7,9,7,3
4,8,7,6,2
5,4,3,2,1
```
All length 7 palindromes
```
1,2,3,4,3,2,1
```

And whatever else you like.

Dependencies:
```
add Folds, StaticArrays, DataStructures, IterTools, Combinatorics
```

Run as such (if you have 16 cores)
```
julia -O3 -t16 Orchestrator.jl square/6 some_az_only_dictionary
```

