qui {

/**********************************************************************************/
/* program masala_merge : Fuzzy match using masalafied levenshtein                */
/*
Meta-algorithm:
1. stata outsheets two files, with an id and a name column.
2. python gets filenames and parameters from command lines, reads two files into two dictionaries
3. python outputs a single file, with id, str1, str2, distance
4. stata reads this file and processes, makes decisions about which matches to keep

See readme.md for sample usage.

*/
/***********************************************************************************/
cap prog drop masala_merge
prog def masala_merge
{
  syntax [varlist] using/, S1(string) OUTfile(string) [FUZZINESS(real 1.0) quietly KEEPUSING(passthru) SORTWORDS] 

  /* require tmp and masala_dir folders to be set */
  if mi("$tmp") | mi("$MASALA_PATH") {
      disp as error "Need to set globals 'tmp' and 'MASALA_PATH' to use this program"
      exit
  }

  /* store sort words parameter */
  if !mi("`sortwords'") {
    local sortwords "-s"
  }
  
  /* define maximum distance for lev.py as 0.35 + 1.25 * (largest acceptable match).
     This is the threshold limit, i.e. if we accept a match at 2.1, we'll reject it
        if there's another match at 0.4 + 2.1*1.25. (this is hardcoded below) */
  local max_dist = 0.40 + 1.25 * 2.1 * `fuzziness'
  
  /* make everything quiet until python gets called -- this output is not helpful */
  qui {

    /* create temporary file to store original dataset */
    tempfile master
    save `master', replace

    /* create a random 5-6 digit number to make the temp files unique */
    local time = real(subinstr(`"`c(current_time)'"', ":", "", .))
    local nonce = floor(`time' * runiform() + 1)
    
    local src1 $tmp/src1_`nonce'.txt
    local src2 $tmp/src2_`nonce'.txt
    local out $tmp/out_`nonce'.txt
    local lev_groups $tmp/lev_groups_`nonce'.dta
    
    preserve
    
    keep `varlist' `s1'
    sort `varlist' `s1'
    
    /* merge two datasets on ids to produce group names */
    merge m:m `varlist' using `using', keepusing(`varlist' `s1')
    
    // generate id groups
    egen g = group(`varlist')
    drop if mi(g)
    
    qui sum g
    local num_groups = r(max)
            
    // save group list
    keep g `varlist'
    duplicates drop
    save "`lev_groups'", replace
    
    /* now prepare group 1 */
    restore
    preserve
    
    keep `varlist' `s1'
  
    /* drop if missing string and store # observations */
    keep if !mi(`s1')
    qui count
    local g1_count = r(N)
    
    /* bring in group identifiers */
    merge m:1 `varlist' using "`lev_groups'", keepusing(g)
  
    /* places with missing ids won't match group */
    drop if _merge == 1
  
    /* only keep matches */
    keep if _merge == 3
    duplicates drop
    
    // outsheet string group 1
    outsheet g `s1' using "`src1'", comma replace nonames
    
    // prepare group2
    di "opening `using'..."
    use `using', clear
    keep `varlist' `s1'

    /* confirm no duplicates on this side */
    tempvar dup
    duplicates tag `varlist' `s1', gen(`dup')
    count if `dup' > 0
    if `r(N)' > 0 {
      display as error "`varlist' `s1' not unique on using side"
      exit 123
    }
    drop `dup'
    
    /* drop if missing string and store # observations */
    keep if !mi(`s1')
    qui count
    local g2_count = r(N)
    
    // merge in group identifiers
    merge m:1 `varlist' using "`lev_groups'", keepusing(g)
    
    /* something wrong if didn't match group ids for any observation */
    drop if _merge == 1
  
    /* only keep matches */
    keep if _merge == 3
    duplicates drop
    
    // outsheet string group 2
    outsheet g `s1' using "`src2'", comma replace nonames
  }
  
  // call python levenshtein program
  di "Matching `g1_count' strings to `g2_count' strings in `num_groups' groups."
  di "Calling lev.py:"

  di `" shell python -u $MASALA_PATH/lev.py -d `max_dist' -1 "`src1'" -2 "`src2'" -o "`out'" `sortwords'"'
  !python               $MASALA_PATH/lev.py -d `max_dist' -1 "`src1'" -2 "`src2'" -o "`out'" `sortwords'

  di "lev.py finished."

  /* quietly process the python output */
  qui {
    /* open output lev dataset */
    /* take care, this generates an error if zero matches */
    capture insheet using "`out'", comma nonames clear
  
    /* if there are zero matches, create an empty outfile and we're done */
    if _rc {
      disp_nice "WARNING: masala_merge: There were no matches. Empty output file will be saved."
      clear
      save `outfile', replace emptyok
      exit
    }
    ren v1 g
    ren v2 `s1'_master
    ren v3 `s1'_using
    ren v4 lev_dist
  
    /* merge group identifiers back in */
    destring g, replace
    merge m:1 g using "`lev_groups'", keepusing(`varlist')
    
    /* _m == 1 would imply that our match list has groups not in the initial set */
    assert _merge != 1
  
    /* _m == 2 are groups with zero matches. drop them */
    drop if _merge == 2
  
    /* count specificity of each match */
    bys g `s1'_master: egen master_matches = count(g)
    bys g `s1'_using: egen using_matches = count(g)
  
    /* count distance to second best match */
  
    /* calculate best match for each var */
    foreach v in master using {
      bys g `s1'_`v': egen `v'_dist_rank = rank(lev_dist), unique
      
      gen tmp = lev_dist if `v'_dist_rank == 1
      bys g `s1'_`v': egen `v'_dist_best = max(tmp)
      drop tmp
      gen tmp = lev_dist if `v'_dist_rank == 2
      bys g `s1'_`v': egen `v'_dist_second = max(tmp)
      drop tmp
      
      drop `v'_dist_rank
    }
    
    drop g _m
  
    /* apply optimal matching rule (based on 1991-2001 pop census confirmed village matches in calibrate_fuzzy.do) */
    /* initialize */
    gen keep_master = 1
    gen keep_using = 1
  
    /* get mean length of matched string */
    gen length = floor(0.5 * (length(`s1'_master) + length(`s1'_using)))
  
    /* 1. drop matches with too high a levenshtein distance (threshold is a function of length) */
    replace keep_master = 0 if lev_dist > (0.9 * `fuzziness') & length <= 4
    replace keep_master = 0 if lev_dist > (1.0 * `fuzziness') & length <= 5
    replace keep_master = 0 if lev_dist > (1.3 * `fuzziness') & length <= 8
    replace keep_master = 0 if lev_dist > (1.4 * `fuzziness') & inrange(length, 9, 14)
    replace keep_master = 0 if lev_dist > (1.8 * `fuzziness') & inrange(length, 15, 17)
    replace keep_master = 0 if lev_dist > (2.1 * `fuzziness')
    
    /* copy these thresholds to keep_using */
    replace keep_using = 0 if keep_master == 0
  
    /* 2. never use a match that is not the best match */
    replace keep_master = 0 if (lev_dist > master_dist_best) & !mi(lev_dist)
    replace keep_using = 0 if (lev_dist > using_dist_best) & !mi(lev_dist)
    
    /* 3. apply best empirical safety margin rule */
    replace keep_master = 0 if (master_dist_second - master_dist_best) < (0.4 + 0.25 * lev_dist)
    replace keep_using = 0 if (using_dist_second - using_dist_best) < (0.4 + 0.25 * lev_dist)
  
    /* save over output file */
    order `varlist' `s1'_master `s1'_using lev_dist keep_master keep_using master_* using_*
    save `outfile', replace
  }
  restore

  /* run masala_review */
  use `outfile', clear
  
  /* if quietly is not specified, call masala_review, which calls masala_process */
  if mi("`quietly'") {
    masala_review `varlist', s1(`s1') master(`master') using(`using')
  }

  /* if quietly is specified, go directly to masala_process */
  else {
    masala_process `varlist', s1(`s1') master(`master') using(`using')
  }
  
  di "Masala merge complete."
  di " Original master file was saved here:   `master'"
  di " Complete set of fuzzy matches is here: `outfile'"
}
end
/* *********** END program masala_merge ***************************************** */

/**********************************************************************************/
/* program masala_lev_dist : Calculate levenshtein distance between two vars */
/*                           uses external python program */
/***********************************************************************************/
cap prog drop masala_lev_dist
prog def masala_lev_dist
{
  syntax varlist(min=2 max=2), GEN(name)
  tokenize `varlist'
  foreach i in _masala_word1 _masala_word2 _masala_dist __masala_merge {
    cap drop `i'
  }

  gen _masala_word1 = `1'
  gen _masala_word2 = `2'
  replace _masala_word1 = lower(trim(_masala_word1))
  replace _masala_word2 = lower(trim(_masala_word2))

  gen _row_number = _n
  
  /* create temporary file for python  */
  outsheet _row_number _masala_word1 _masala_word2 using $tmp/masala_in.csv, comma replace nonames

  /* call external python program */
  di "Calling lev.py..."
  shell python $MASALA_PATH/lev.py -1 $tmp/masala_in.csv -o $tmp/masala_out.csv

  /* convert created file to stata format */
  preserve
  insheet using $tmp/masala_out.csv, clear names
  save $tmp/masala_lev_dist, replace
  restore

  /* merge result with new dataset */
  merge 1:1 _row_number using $tmp/masala_lev_dist.dta, gen(__masala_merge) keepusing(_masala_dist)

  /* clean up */
  destring _masala_dist, replace
  ren _masala_dist `gen'
  drop _masala_word1 _masala_word2 _row_number
  
  assert __masala_merge == 3
  drop __masala_merge
}
end
/* *********** END program masala_lev_dist ***************************************** */

/*****************************************************************************************************/
/* program fix_spelling : Fixes spelling in a string variable based on a supplied master key.        */
/*                Runs a fuzzy masala-merge of data-list to master-list within group if specified    */

/* Process: take data in data_list, merge it to master list
            keep if _m == 2
            fuzzy merge data-list to master-list, maybe within some group.
            if a single match without competition, then replace data list with the data in master list
            return new version of data list                                                          */
/* Must specify gen() or replace                                                                     */

/* Targetfield and targetgroup() allow varnames in the master set to differ from varnmes in the key  */
/* also allow running of program for pc01_district_name if a variable pc01_district_name already exists
   in set
targetfield:  refers to variable name in key named according to pc`year'_`place'_name
targetgroup:  refers to variable names of group names on which to group spelling replacements, i.e.
              i.e. group is pc01_state_name if fix_spelling pc01_district_name                       */
/* Example syntax:

fix_spelling district_name_fm, targetfield(pc01_district_name) group(state_name_fm) targetgroup(pc01_state_name) ///
  src($keys/pc01_district_key) gen(new_district_name)

where target is universal naming syntax in keys, and state_name_fm / district_name_fm are names in set */
/*******************************************************************************************************/
cap prog drop fix_spelling
prog def fix_spelling
{
  syntax varname(min=1 max=1), SRCfile(string) [GEN(name) GROUP(varlist) TARGETfield(string) TARGETGROUP(string) keepall replace FUZZINESS(real 2)]

  /* put variable fix spelling in `varname' because we need arg "`1'" open */
  tokenize `varlist'
  local varname `1'
  
  /* need to specify either generate or replace, not both */
  if (mi("`gen'") & mi("`replace'")) | (!mi("`gen'") & !mi("`replace'")) {
    display as error "fix_spelling: Need to specify either generate or replace, not both"
    exit 1
  }

  /* if replace is set, create a temp var to be generated */
  if !mi("`replace'") {
    tempvar gen
  }
  
  /* if group is empty, need to create a master group that is the entire file */
  if mi("`group'") {
    tempvar __GROUP
    gen `__GROUP' = 1
    local nogroup = 1
    local group "`__GROUP'"
  }
  
  /* for now, assume we have a source file */

  /* create the master list */
  qui {
    preserve

    /* open the synonym list: master location key */
    use "`srcfile'", clear

    /* if TARGETFIELD was specified, rename the target to the variable we want to match */
    if !mi("`targetfield'") {
      ren `targetfield' `varname'
    }

    /* if TARGETGROUP was specified, rename group variables to match key used to spellcheck */
    if !mi("`targetgroup'") {
      
       /* loop through each element in group and targetgroup to rename each variable from targetgroup as named in master set */
       /* group -> tokens `1' -> `n', while we loop over target group. Need to do this because can't tokenize two strings */
       tokenize `group'
       local i = 1
       foreach targetgroup_var in `targetgroup' {
         cap confirm variable `targetgroup_var'
         if _rc {
           disp as error "ERROR fix_spelling: key file missing targetgroup var `targetgroup_var'"
           exit 123
         }
         ren `targetgroup_var' ``i'' 
         local i = `i' + 1
       }         
     }
    
    /* define a group if none exists */
    if !mi("`nogroup'") gen `__GROUP' = 1
    
    /* save clean, unique synonym list as stata file */
    tempfile SPELLING_MASTER_LIST
    keep `group' `varname'
    duplicates drop
    sort `group' `varname'
    save `SPELLING_MASTER_LIST', replace
    restore
    
    /* create a list of unmatched names */
    preserve
  
    keep `varname' `group'
    duplicates drop
  
    /* get rid of exact matches - these will work well */
    merge 1:1 `group' `varname' using `SPELLING_MASTER_LIST', gen(_merge1)
    keep if _merge1 == 1
  
    /* if nothing left, then the original list is fine and we're done */
    qui count
    if r(N) == 0 {
      restore
      noi di "100% of names matched exactly. No fuzzy matching necessary"
      
      /* pass variable back into original or specified gen variable */
      gen `gen' = `varname'
      
      /* drop group var from main set if no group was specified */
      if !mi("`nogroup'") drop `__GROUP'
      exit
    }

    /* set tempfile to outfile spelling errors */
    tempfile spelling_errors

    /* otherwise, go to the fuzzy merge */
    masala_merge `group' using `SPELLING_MASTER_LIST', s1(`varname') outfile(`spelling_errors') fuzziness(`fuzziness')

    /* set tempfile for spelling corrections */
    tempfile SPELLING_CORRECTIONS
    
    /* review masala merge results */
    use `spelling_errors', clear
  
    /* exit if no matches */
    count
    if `r(N)' == 0 exit
  
    /* keep best match for everything in badly-spelled set */
    keep if keep_master == 1
    keep `group' `varname'_master `varname'_using lev_dist
  
    /* fix names and merge back to the original dataset */
    ren `varname'_master `varname'
    ren `varname'_using `gen'
    ren lev_dist `gen'_dist
    save `SPELLING_CORRECTIONS', replace
    restore
  
    /* tag exact matches (this merge only adds _merge_exact) */
    merge m:1 `group' `varname' using `SPELLING_MASTER_LIST', gen(_merge_exact)
    drop if _merge_exact == 2
    
    /* then get fuzzy matches */
    merge m:1 `group' `varname' using `SPELLING_CORRECTIONS', gen(_merge_fuzzy)
    assert _merge_fuzzy != 2
  
    /* if we have an exact match, shouldn't have a fuzzy match */
    assert _merge_fuzzy == 1 if _merge_exact == 3
    
    /* add exact matches */
    replace `gen' = `varname' if _merge_exact == 3
    replace `gen'_dist = 0 if _merge_exact == 3
    drop _merge_exact _merge_fuzzy
  
    /* if keepall specified, get places that didn't match */
    if !mi("`keepall'") {
  
      /* merge the spell-checked data back to the master list within the group */
      ren `varname' `varname'_SP
      ren `gen' `varname'
      merge m:1 `group' `varname' using `SPELLING_MASTER_LIST', nogen keepusing(`varname')
      ren `varname' `gen'
      ren `varname'_SP `varname'
    }

    /* if group was not specified (or place is state so there is no group), drop group var */
    if !mi("`nogroup'") drop `__GROUP'
  }
  /* if replace was specified */
  if !mi("`replace'") {

    /* show replacements made */
    tempvar tag
    qui egen `tag' = tag(`varname') if !mi(`gen') & `gen' != `varname'
    disp_nice "Spelling fixes and Masala-Levenshtein distances:"
    list `varname' `gen' `gen'_dist if `tag'
    
    /* replace original var, show what was done, and drop the distance */
    qui replace `varname' = `gen' if !mi(`gen')
    drop `gen' `gen'_dist
  }

}
end
/* *********** END program fix_spelling ***************************************** */

/**********************************************************************************/
/* program masala_review : Reviews masala_merge results and calls masala_process  */
/***********************************************************************************/
cap prog drop masala_review
prog def masala_review
{
  syntax varlist, s1(string) master(string) using(string) [keepusing(passthru)]

  /* ensure a masala merge output file is open */
  cap confirm var keep_master
  if _rc {
    di "You must open the masala_merge output file before running this program."
  }
  
  /* count and report matches that are exact, but with alternatives */
  /* these are places where keep_master == 0 & lev_dist == 0 */
  qui bys `s1'_master: egen _min_dist = min(lev_dist)
  qui bys `s1'_master: egen _max_keep = max(keep_master)

  qui count if _max_keep == 0 & _min_dist == 0
  if `r(N)' > 0 {
    di "+-------------------------------------" _n "| These are exact matches, where alternate good matches exist." _n ///
      "| keep_master is 0, but masala_process() will keep the observations with lev_dist == 0." _n ///
        "+-------------------------------------" 
    list `varlist' `s1'* lev_dist if _max_keep == 0 & _min_dist == 0
  }
  qui drop _max_keep _min_dist

  /* visually review places with high lev_dist that script kept -- they look good. */
  qui count if keep_master == 1 & lev_dist > 1
  if `r(N)' > 1 {
    disp_nice "These are high cost matches, with no good alternatives. keep_master is 1."
    list `varlist' `s1'* lev_dist if keep_master == 1 & lev_dist > 1
  }

  /* run masala_process, and then show the unmatched places */
  masala_process `varlist', s1(`s1') master(`master') using(`using') `keepusing'

  /* tag each name so it doesn't appear more than once */
  qui egen _ntag = tag(`varlist' `s1')

  /* list unmatched places in a nice order */
  qui gen _matched = _masala_merge == 3
  gsort _matched -_ntag `varlist' `s1'

  /* ensure we don't trigger obs. nos. out of range in final list, by counting observations */
  qui count
  if `r(N)' < 200 {
    local limit `r(N)'
  }
  else {
    local limit 200
  }

  /* list unmatched places */
  qui count if _masala_merge < 3 & _ntag in 1/`limit'
  if `r(N)' {
    disp_nice "This is a sorted list of some places that did not match. Review for ideas on how to improve"
    list `varlist' `s1' _masala_merge if _masala_merge < 3 & _ntag in 1/`limit', sepby(`varlist')
  }

  drop _ntag _matched
}
end
/* *********** END program masala_review ***************************************** */

/**********************************************************************************/
/* program masala_process : Rejoins the initial files in a masala_merge           */
/**********************************************************************************/
cap prog drop masala_process
prog def masala_process
{
  syntax varlist, s1(string) master(string) using(string) [keepusing(passthru)]

  qui {
    /* override keep_master if lev_dist is zero. */
    replace keep_master = 1 if lev_dist == 0
    
    /* keep highly reliable matches only */
    keep if keep_master == 1
  
    /* drop all masala merge's variables */
    keep `varlist' `s1'*

    /* bring back master dataset */
    gen `s1' = `s1'_master
    merge 1:m `varlist' `s1' using `master', gen(_masala_master)

    /* fill in master fuzzy-string from unmatched data on master side */
    replace `s1'_master = `s1' if mi(`s1'_master)
    drop `s1'
    
    /* bring back using dataset */
    gen `s1' = `s1'_using
  
    merge m:1 `varlist' `s1' using `using', `keepusing' gen(_masala_using)

    /* fill in using fuzzy-string from unmatched data on using side  */
    replace `s1'_using = `s1' if mi(`s1'_using)
    drop `s1'
    
    /* set `s1' to the master value */
    ren `s1'_master `s1'
    
    /* fill in using values when _m == 2 */
    replace `s1' = `s1'_using if mi(`s1')
  }

  /* Assertion: if we couldn't match back to the using, it must be unmatched from the master side */
  assert _masala_master == 2 if _masala_using == 1

  /* show merge result */
  disp_nice "Results of masala_merge (counting unique strings only): "
  
  /* tag each name so it doesn't appear more than once */
  qui egen ntag = tag(`varlist' `s1')

  /* create a standard merge output variable */
  qui gen _masala_merge = 1 if _masala_master == 2
  qui replace _masala_merge = 2 if _masala_using == 2
  qui replace _masala_merge = 3 if _masala_using == 3 & _masala_master == 3
  drop _masala_master _masala_using
  label values _masala_merge _merge

  /* show results */
  table _masala_merge if ntag
  qui drop ntag
}
end
/* *********** END program masala_process ***************************************** */

/**********************************************************************************/
/* program review_merge : call right after a merge to review potential matches    */
/***********************************************************************************/
cap prog drop review_merge
prog def review_merge
{
  syntax varlist [if/], [merge(string) list(varlist) SHOWall]

  tempvar tag
  
  if mi("`merge'") {
    local merge _merge
  }

  /* `showall' parameter determines whether we limit to _merge < 3 or show all results */
  if !mi("`showall'") {
    local show 1
  }
  else {
    local show `merge' < 3
  }

  /* need something in `if' if nothing passed in */
  if mi("`if'") {
    local if 1
  }
  
  /* separate list with sepby() if more than one variable in varlist */
  tokenize "`varlist'"
  if !mi("`2'") {
    local sepby sepby(`1')
  }

  /* only show each posisble match once */
  egen `tag' = tag(`varlist')
  
  sort `varlist' `list'
  list `varlist' `list' `merge' if `show' & `if' & `tag', `sepby'
  drop `tag'
}
end
/* *********** END program review_merge ***************************************** */

/***********************************************************************************/
/* program create_merge_fragments : Call right after a merge to create separate
                                    files of the unmatched pieces                  */
/***********************************************************************************/
cap prog drop create_merge_fragments
prog def create_merge_fragments
{

  /* idea is we had:
 file1: merge_vars a b c
 file2: merge_vars d e f

we want to create file1 and file2 leftovers. hard part is getting the variables right.

syntax option:
- call with the completed merge file open, pass original master() and using() files back in.

  */
  syntax anything, master(string) using(string) [merge(string) suffix(string)]

  /* set default values for merge and suffix locals */
  if mi("`merge'") local merge _merge
  if mi("`suffix'") local suffix unmatched
  
  /* hack work to get m:1, 1:1, or 1:m */
  local merge_type = substr("`anything'", 1, 3)

  if !inlist("`merge_type'", "1:1", "m:1", "1:m") {
    di "Must specify 1:1, m:1 or 1:m as in merge syntax"
    barf
  }
  local master_type = substr("`merge_type'", 1, 1)
  local using_type = substr("`merge_type'", 3, 1)
  
  local varlist = substr("`anything'", 4, .)
  
  /* confirm varlist is a varlist */
  confirm var `varlist'

  /* we want to leave this unaltered, so wrap everything in a preserve / restore */
  preserve

  /* keep only the matches and drop _merge */
  keep if `merge' == 3
  drop `merge'
  
  /* save the file with the matches. all we need is the varlist */
  keep `varlist'

  /* we only need one copy of each match, this allows everything below to be m:1 */
  duplicates drop
  tempfile merge3
  save `merge3', replace
  
  /* create master only file */
  use `master', clear

  /* merge it to list of matches */
  /* 1:m is con->village. merged file will have many repeated cons */
  /* m:1 is village->con. merged file will have each village once */
  /* 1:1 obviously has each side once */
  merge 1:`using_type' `varlist' using `merge3'

  /* now we want to keep only the non-matches */
  keep if `merge' == 1

  /* there should not be any using, if we just ran this merge */
  assert `merge' != 2

  /* drop _merge and save fragment file */
  drop `merge'
  save `master'_`suffix', replace

  /* repeat process for using side */
  use `using', clear

  /* merge it to list of matches */
  /* 1:m is con->village. merged file will have each village once */
  /* m:1 is village->con. merged file will have cons repeated */
  /* 1:1 obviously has each side once */
  merge `master_type':1 `varlist' using `merge3'

  /* now we want to keep only the non-matches */
  keep if `merge' == 1

  /* there should not be any using, if we just ran this merge */
  assert `merge' != 2

  /* drop _merge and save fragment file */
  drop `merge'
  save `using'_`suffix', replace
  
  restore

  /* report what happened */
  di "Created files with merge fragments:"
  di "  master: `master'_`suffix'"
  di "  using: `using'_`suffix'"
}
end
/* *********** END program create_merge_fragments ***************************************** */

/*****************************************************************************************************/
/* program synonym_fix : Replace strings for names that have multiple variants but one main version, */
/*                       using a externally-supplied dictionary                                      */
/*                       i.e. uttaranchal -> uttarakhand                                             */
  /* external file has the following columns:
      "master"             : the target name
      "pc01_district_name" : variable that we are replacing
      "`group'"            : a list of variables that define the context in which to match the name
                             named according to standard pc`year'_`place'_name
   external file must be in CSV format.
*/
/* options targetfield and targetgroup allow variables in dataset to differ from variables in .csv
targetfield:  refers to column name in external csv file named according to pc`year'_`place'_name
targetgroup:  refers to the columns after the main name column on which to group replacements, i.e.
              synonym_fix pc01_district_name has group pc01_state_name
*/

/*****************************************************************************************************/
cap prog drop synonym_fix
prog def synonym_fix
{
  syntax varname, SYNfile(string) [GENerate(name) replace group(varlist) TARGETfield(string) TARGETGROUP(string) INSHEET] 
  
  /* put variable to replace using synonym_fix in `varname' */
  tokenize `varlist'
  local varname `1'
  
  /* require generate or replace */
  if (mi("`generate'") + mi("`replace'")) != 1 {
    display as error "synonym_fix: generate or replace must be specified"
    exit 1
  }
  
  /* if no generate specified, make replacements to passed in variable */
  if mi("`generate'") {
    local name = "`varname'"
  }

  /* if generate specified, copy the variable into the new slot and then change it gradually */
  else {
    gen `generate' = `varname'
    local name = "`generate'"
  }
  
  qui {
    
    /* verify 'master' variable is not in use in order to make replacements */
    /* master refers to first column in csv that contains the correct string replacement value */
      cap confirm variable master
      if !_rc {
        display as error "'master' variable already exist. synonym_fix() needs this varname to be free."
        exit 123
      }
      
      /* insheet external csv from ~/iecmerge/pc[]/place_synonyms/  */
      preserve

      /* if insheet is specified, read the synfile by using insheet using */
      if !mi("`insheet'") {
        insheet using "`synfile'", clear delim(",") name
      }

      /* if insheet is not specified, read the synfile by using import delimited  */
      else {
        import delimited "`synfile'", clear delim(",") varn(1)
      }

      /* targetfield: renaming the target name variable from the replacement file to match the master dataset */

      /* if a target variable field was specified, rename the insheeted target field variable to match the varname in master dataset */
      if !mi("`targetfield'") {
        cap confirm variable `targetfield'
        if _rc {
          disp as error "ERROR synonym_fix: input file missing target var `targetfield'"
          exit 123
        }
        ren `targetfield' `varname'
      }

      /* otherwise, if these names are the same, do a clean check to confirm synfix has the right field */
      else {
        cap confirm variable `varname'
        if _rc {
          disp as error "ERROR synonym_fix: input file missing variable `varname'"
          exit 123
        }
      }
      
     /* targetgroup: renaming group variables from replacement file to match master dataset */

     /* if a target group was specified, rename the target group names from csv to match master group varnames passed in with group() */
     /* if targetgroup is specified it is implicit that group has been specified */
      if !mi("`targetgroup'") {
        
        /* loop through each element in group and targetgroup to rename each variable from targetgroup as named in master set */
       /* group -> tokens `1' -> `n', while we loop over target group. Need to do this because can't tokenize two strings */
        tokenize `group'
        local i = 1
        foreach targetgroup_var in `targetgroup' {
          cap confirm variable `targetgroup_var'
          if _rc {
            disp as error "ERROR synonym_fix: input file missing targetgroup var `targetgroup_var'"
            exit 123
          }
          ren `targetgroup_var' ``i'' 
          local i = `i' + 1
        }         
      }

      /* assert no duplicates on replacement value */
      cap bys `varname' `group': assert _N == 1
      if _rc {
        display as error ""
        display as error "ERROR synoynm_fix: Synonym file must be unique on value to be replaced (column 2) and optional merge variables."
        noi di "Target group was: `group'"
        noi duplicates list `varname' `group'
        exit 123
      }
      
      /* write the dictionary of replacements */
      tempfile dict
    
      /* the csv should already be formatted lower and trimmed but perform replace in case it's not */
      replace master = trim(lower(master))
      foreach item in `group' {
        replace `item' = trim(lower(`item'))
      }
    
      /* drop empty rows in case stata insheeted blank rows from the csv */      
      drop if mi(master)
      save "`dict'", replace
      restore
      
      /* prepare group vars for merge */
      if !mi("`group'") {
        foreach v of varlist(`group') {
          tostring `v', replace force
          replace `v' = trim(lower(`v'))
        }
      }
      
      /* merge passed in variable to the synonym column */
      merge m:1 `group' `varname' using `dict', keepusing(master)
      drop if _merge == 2
    }    

    /* make the actual string replacement - put out of quietly block so user can see what happened */
    replace `name' = master if !mi(master)

    drop master _merge
  }
end
/* *********** END program synonym_fix ***************************************** */

/********************************************************/
/* program masala_merge2 : Placeholder for masala_merge */
/********************************************************/
cap prog drop masala_merge2
prog def masala_merge2
{
  syntax [varlist] using, S1(passthru) OUTfile(passthru) [FUZZINESS(passthru) quietly(passthru) KEEPUSING(passthru) SORTWORDS(passthru)] 
  masala_merge `varlist' `using', `s1' `outfile' `fuzziness' `quietly' `keepusing' `sortwords' `quietly'
}
end
/* *********** END program masala_merge2 ***************************************** */


/**********************************************************************************/
  /* program disp_nice : Insert a nice title in stata window */
  /***********************************************************************************/
  cap prog drop disp_nice
  prog def disp_nice
    {
      di _n "+--------------------------------------------------------------------------------------" _n "| `1'" _n  "+--------------------------------------------------------------------------------------"
    }
  end
  /* *********** END program disp_nice ***************************************** */

}

