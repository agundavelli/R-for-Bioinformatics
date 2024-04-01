[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-24ddc0f5d75046c5622901739e7c5dd533143b0c8e959d652212380cedb1ea36.svg)](https://classroom.github.com/a/Qkqy67w_)
# Assignment Differential Expression Part 2

> Volcanoes don't really look like that, do they?

## Problem Statement
Because R has been utilized in bioinformatics communities for years at this point, there are many different approaches to do one task. While in the last assignment we just used DESeq2, we will be performing differential expression analysis with two additional packages (edgeR and Limma) and comparing the results we find.

## Learning Objectives
1. The basic operations of `ggplot` and "The Grammar of Graphics".
2. Using differential gene expression analysis packages to generate data for plotting.
3. Creating simple plots and **not** using the default colors and themes.
4. Combining multiple plots into one image using `facet_wrap()`

## Skill List
- A tempered heart and mind that understands a plot may not ever look _exactly_ how you want it to.
- An intermediate understanding of R's most popular plotting package, `ggplot`.
- Further understanding of differential expression analysis and its plotting.
- A sense of superiority whenever you see a publication use the base R plotting package for figures.

## Instructions
Complete `main.R` in a way that satisfies all of the tests in `test_main.R`. Follow the instructions in the function descriptions for `main.R` and read here for details on the tests your code should satisfy. You will use these functions to complete figures in `report.Rmd`.

You can use `testthat:test_file('test_main.R')` in your R console to run all tests in `test_main.R`. Do not modify the tests. While it is difficult to write tests that cover every solution to a problem, our tests are written to ensure the outputs of a function are aligned with what an analysis should perform. If you feel a test is working incorrectly, let a TA know.

## Function Details

### 1. `load_n_trim()`
Once again we have some data to load into R in order to manipulate. We'll be loading a counts file, which is a matrix of genes (rows) by samples (columns). Each cell is the number of that gene found in that sample. This is a good opportunity to use a data frame or a tibble, but because of package restrictions we will be sticking to data frames for this assignment. We're using a couple different packages to process our data, and they input/output _exclusively_ in data frames, so we will be consistent and save ourselves some trouble with conversions.  

We want to return a data frame that is essentially identical to our input file (always a good idea to use the command line to see what you're working with). Ensure the gene names are stored as `row.names` and _not_ as a separate column, data frames have a distinction between the two. This function should also reduce the columns to just those of interest (P0 and Adult).

***Hints***
Row names can be changed using the `row.names()` function. The new row names must be a character vector of exactly the same length.
``` 
df <- data.frame(a = c(0, 0, 0), b = c(1, 3, 5))
row.names(df) <- c("Row1", "Row2", "Row3")
print(df)
     a b
Row1 0 1
Row2 0 3
Row3 0 5
```


#### Tests
```
describe("Data Loading and Reduction", {
  
  test_df <- load_n_trim(csv)
  
  it("should have the correct column names", {
    expect_equal(names(test_df), c("vP0_1", "vP0_2", "vAd_1", "vAd_2"))
  })
  
  it("should have the correct number of rows and columns", {
    expect_equal(dim(test_df), c(55416, 4))
  })
  
  it("should be a data.frame", {
    expect_equal(class(test_df), "data.frame")
  })
})
```

The tests for this assignment are relatively straightforward, mostly ensuring that data is loaded and ends up in the right shape and format. The first test ensures the columns were filtered, the second ensures the size is correct, and the final ensures it is a data frame and _not_ a tibble. No tibbles. Tibble-less?

### 2. `run_deseq()`
This is the first of three functions that technically do identical things (but in different ways). One of the most popular Bioconductor packages, DESeq2 has a number of options available for differential expression analysis. We will load in the counts data, select our variables of interest, and use `DESeq()` to process them. Links are included in the function description to DESeq2 documentation, you will find answers to most questions there, especially since those were the documents used to _write_ this assignment. This function will return the results of the analysis.

#### Tests
```
describe("DESeq2 Functionality", {
  
  load("data/mock_counts_df.RData") # loads the counts_df object into env
  coldata <- data.frame(condition = rep(c("day4", "day7"), each=2),
                        type="paired-end")
  row.names(coldata) <- c("vP4_1", "vP4_2", "vP7_1", "vP7_2")
  deseq <- run_deseq(counts_df, coldata, 10, "condition_day7_vs_day4")
  
  it("should return a DESeqResults object", {
    expect_equal(class(deseq)[1], "DESeqResults")
  })
  
  it("should have the correct dimensions", {
    expect_equal(dim(deseq), c(19127, 6))
  })
  
  it("should have pvalue and padj columns", {
    expect_equal(c("pvalue", "padj") %in% names(deseq), c(TRUE, TRUE))
  })
})
```

Less concerning for you all but writing tests for larger data sets and more complicated functions like this is...tricky. The problem is it can't be written like a script, we can't use the results from `load_n_trim()` to test _this_ function, because what if one fails independently of the other? So in this case, we create an RData object of correctly loaded data that works _even if_ `load_n_trim()` is absolutely broken. Also, in this case, we are using different sets to test this (days 4 and 7 instead of 0 and adult) so you can't just load the test data and say you wrote the code! That would be disingenuous of you, and being disingenuous is for after you graduate.

This test loads the sample data (a data frame of counts) and creates the `coldata` object according to DESeq2 specifications. We wrap the `run_deseq()` in `expect_warning()` because DESeq throws a warning when it performs an expected behavior: converting strings to factors (this is bad programming! this is not what warnings are for, the package is functioning as expected). We then do a similar series of tests: ensure the dimensions are correct, ensure the results object are from DESeq, and ensure that the column names we need are present.

### 3. `run_edger()`
As mentioned in the previous function, this is the EdgeR implementation of running the differential expression analysis based on our counts file. Again, this function should follow roughly the portion of the documentation mentioned in the function description. EdgeR has a number of additional plotting functions (such as `plotMDS()`) you may like to include. These are optional but pretty fun to see if your data is working out as expected (similar days clustered close together and so on). You will return the _entire_ results object, no filtering necessary yet.

#### Tests
```
describe("edgeR Functionality", {
  
  load("data/mock_counts_df.RData")
  group <- factor(rep(c(1,2), each=2))
  edger_res <- run_edger(counts_df, group)
  
  it("should have the correct column names", {
    expect_equal(names(edger_res), c("logFC", "logCPM", "PValue"))
  })
  
  it("should have the correct dimensions", {
    expect_equal(dim(edger_res), c(15026, 3))
  })
})
```

Similar to the DESeq2 tests, we aren't necessarily looking for a carbon copy of the results. While the process _should_ be entirely deterministic, it is better to allow some wiggle room and merely look that the columns we're interested in are present and that the data has the right shape. Note that EdgeR doesn't perform a BH correction (the p-adjusted column is non-existant). It's either an option that's off by default or just not a feature! Not a problem, we have the p-values so we can do it ourselves (see the markdown).

### 4. `run_limma()` and `run_voom()`
Our final in this set of three, we're testing two related functions in the Limma package. While you can create an analysis that _only_ uses Limma, in our case we can apply the **voom** workflow to refine our data a little more. More info on how these work in the cited documentation. A little technique that can be useful for optional workflows like this is including a boolean flag in the parameters of the function. If I am running a function `add_numbers(1, -4, abs=TRUE)`, the parameter `abs` might change the functions operation so that it takes the absolute value of numbers before adding them together. I could create this behavior like so:
```
add_numbers <- function(x, y, abs=FALSE) {
  if (abs) { # code enters here if abs is TRUE
    return(abs(x) + abs(y)) # return means function exits, doesn't go to "else {}"
  } else { # code enters here if abs was false
    return(x + y)
  }
}
```
This is a very useful technique for re-using functions, meaning I don't need to write two functions that do similar things when I can change their behavior with that boolean parameter.

#### Tests
```
describe("Limma + Voom Integration", {
  
  load("data/mock_counts_df.RData")
  group <- factor(rep(c(1,2), each=2))
  design <- data.frame(day4=1, day4vsday7=c(0, 0, 1, 1))
  row.names(design) <- c("vP4_1", "vP4_2", "vP7_1", "vP7_2")
  voom_res <- run_limma(counts_df, design, group)
  
  it("should have the correct dimensions", {
    expect_equal(dim(voom_res), c(15026, 6))
  })
  
  it("should have the correct column names", {
    expect_equal(names(voom_res), c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"))
  })
})
```

This ultimately functions similarly to the first two series of tests: load test data, apply the function to it, make sure the results are the right size, shape, and type. This includes the voom procedure, so using _only_ limma results will yield failures.  

***Hints***
Note that Limma uses EdgeR to input data, so your first couple steps will look pretty familiar.


### 5. `combine_pval()`
The final three functions concern modifying our differential expression results and creating plots of these combined results. While it would be trivial to make three plots using the same code or function, have R generate one consistent plot after changes in your data or method is extremely useful.

We use `combine_pval()` to coerce our data into long format, which is more applicable for facet wrapping. `facet_wrap()` is a `ggplot2` function that uses a property of the data to create separate plots in one image. In our case we will use the `package` property (DESeq2, edgeR, or Limma) to differentiate our data. There are a couple ways to reorganize data in this way, but I would suggest using `tidyr::gather()`. Our `key` will be the package and the `value` will be the p-value.


### 6. `create_facets()`
Now that we have examined the p-values, we can look at the adjusted p-values. We'll be constructing volcano plots, which compare the log<sub>2</sub> fold-change and the adjusted p-value, so the `create_facets()` function will be somewhat more complicated. Rather than use `gather()` again, it may be easier to create three tibbles or data frames of the variables of interest and combine them with `rbind()`.


### 7. `theme_plot()`
Finally, we have all of our data in the right shape and we can do some plotting. There is a lot to learn about ggplot, and I mention a few resources in the function description, but don't be afraid to experiment! Try different colors, themes, and styles. The easiest part is getting your data on the plot, the hard part is making it look attractive.

Once all your tests are passing, take a look at the report and see how everything can work together!
