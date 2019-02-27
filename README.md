# masala-merge
Stata/Python code for fuzzy matching of latin script location names in Hindi.

# Stata usage (masala_merge)

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

# Stata usage (fix_spelling)

`fix_spelling` is useful if you want to rationalize your names to a
master list. For example, suppose you have a dataset with district
names, you have a master list of district names (with state
identifiers), and you want to modify your current district names to
match the master key. You could run the following command:

```
fix_spelling district_name, src(master_district_key.dta) group(state_id) replace
```

Additional options:

- `gen(varname)` can be used instead of replace
- `targetfield()` and `targetgroup()` can be used if the group or merge
  variable have different names in the master dataset.
- If `keepall` is specified, your dataset will add rows from the
target list that didn't match anything in your data

Example:

```
. fix_spelling pc01_district_name, src($keys/pc01_district_key) replace group(pc01_state_name) 

[...]

+--------------------------------------------------------------------------------------
| Spelling fixes and levenshtein distances:
+--------------------------------------------------------------------------------------

      +-------------------------------------------+
      | pc01_distri~e         __000000   __0000~t |
      |-------------------------------------------|
  80. |   karimanagar       karimnagar         .2 |
 155. |  mahabubnagar      mahbubnagar         .2 |
 422. |         buxor            buxar        .45 |
 462. |     jahanabad        jehanabad        .45 |
 480. |     khagari a         khagaria        .01 |
      |-------------------------------------------|
 544. |        purnea           purnia        .45 |
 624. |     ahmedabad        ahmadabad        .45 |
 700. |   banaskantha     banas kantha        .01 |
 757. |         dahod            dohad         .8 |
 888. |    panchmahal     panch mahals       1.01 |
      |-------------------------------------------|
 932. |   sabarkantha     sabar kantha        .01 |
 991. |       vadodra         vadodara         .2 |
1490. |         angul           anugul          1 |
1546. |         boudh            baudh        .45 |
1569. |       deogarh         debagarh        1.2 |
      |-------------------------------------------|
1609. | jagatsinghpur   jagatsinghapur         .2 |
1617. |        jajpur          jajapur         .2 |
1674. |        khurda          khordha        .65 |
1722. |   nabarangpur     nabarangapur         .2 |
1922. |    puducherry      pondicherry       1.35 |
      +-------------------------------------------+
```
