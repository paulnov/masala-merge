# masala-merge
Stata/Python code for fuzzy matching of latin script location names in Hindi.

# Stata usage

You should use `masala_merge2`, not `masala_merge`. The latter is included
for backward compatibility with our codebase.

Suppose you already have datasets that align on state, district and
subdistrict codes, and you want to fuzzy match on village names:

```
masala_merge2 state_id district_id subdistrict_id using match_data.dta, S1(village_name) OUTfile($tmp/foo)
```

The outfile is a temporary file that will be generated with the full
set of matches. `masala_merge2` will create a dataset that
automatically picks the matches that it thinks are best and rejects
the matches that are too ambiguous.

Some Additional options:
```
fuzziness(real): Default is 1. Higher numbers generate merges with more
           tolerance for differences. But even 0 will allow some
	   fuzziness.

quietly: Suppress output

keepusing(string): Works like keepusing() in standard merge command

sortwords: Will sort words before running merge
```

