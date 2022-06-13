# What do quality scores mean?

## Phred Quality Scores

We can calculate the *phred quality score* from the probability of
sequencing error (i.e. the base call is wrong) using:

Alternatively, we can rearrange to calculate the probability of error
from the *phred quality score* using:

Where *Q* is the *phred quality score* and *p* is the probability of
error (i.e. the probability that the base call is wrong)

<!-- #endregion -->

``` r
library(tibble)
library(knitr)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
error_prob = function(quality){
  return (10**(quality/-10.0))
}

tibble(`Phred Score`=0:41,
       `Prob of Error`=error_prob(`Phred Score`)) ->
  phred_quality


phred_quality %>%
  kable
```

| Phred Score | Prob of Error |
|------------:|--------------:|
|           0 |     1.0000000 |
|           1 |     0.7943282 |
|           2 |     0.6309573 |
|           3 |     0.5011872 |
|           4 |     0.3981072 |
|           5 |     0.3162278 |
|           6 |     0.2511886 |
|           7 |     0.1995262 |
|           8 |     0.1584893 |
|           9 |     0.1258925 |
|          10 |     0.1000000 |
|          11 |     0.0794328 |
|          12 |     0.0630957 |
|          13 |     0.0501187 |
|          14 |     0.0398107 |
|          15 |     0.0316228 |
|          16 |     0.0251189 |
|          17 |     0.0199526 |
|          18 |     0.0158489 |
|          19 |     0.0125893 |
|          20 |     0.0100000 |
|          21 |     0.0079433 |
|          22 |     0.0063096 |
|          23 |     0.0050119 |
|          24 |     0.0039811 |
|          25 |     0.0031623 |
|          26 |     0.0025119 |
|          27 |     0.0019953 |
|          28 |     0.0015849 |
|          29 |     0.0012589 |
|          30 |     0.0010000 |
|          31 |     0.0007943 |
|          32 |     0.0006310 |
|          33 |     0.0005012 |
|          34 |     0.0003981 |
|          35 |     0.0003162 |
|          36 |     0.0002512 |
|          37 |     0.0001995 |
|          38 |     0.0001585 |
|          39 |     0.0001259 |
|          40 |     0.0001000 |
|          41 |     0.0000794 |

## Ascii Codes

In FASTQ files, phred scores are represented using characters. Each
character on the keyboard can be represented by a number, called an
ascii code.

``` r
tibble(`ASCII #`=33:90) %>%
  rowwise %>%
  mutate(Character=intToUtf8(`ASCII #`)) %>%
  kable
```

| ASCII \# | Character |
|---------:|:----------|
|       33 | !         |
|       34 | "         |
|       35 | \#        |
|       36 | $         |
|       37 | %         |
|       38 | &         |
|       39 | ’         |
|       40 | (         |
|       41 | )         |
|       42 | \*        |
|       43 | \+        |
|       44 | ,         |
|       45 | \-        |
|       46 | .         |
|       47 | /         |
|       48 | 0         |
|       49 | 1         |
|       50 | 2         |
|       51 | 3         |
|       52 | 4         |
|       53 | 5         |
|       54 | 6         |
|       55 | 7         |
|       56 | 8         |
|       57 | 9         |
|       58 | :         |
|       59 | ;         |
|       60 | \<        |
|       61 | =         |
|       62 | \>        |
|       63 | ?         |
|       64 | @         |
|       65 | A         |
|       66 | B         |
|       67 | C         |
|       68 | D         |
|       69 | E         |
|       70 | F         |
|       71 | G         |
|       72 | H         |
|       73 | I         |
|       74 | J         |
|       75 | K         |
|       76 | L         |
|       77 | M         |
|       78 | N         |
|       79 | O         |
|       80 | P         |
|       81 | Q         |
|       82 | R         |
|       83 | S         |
|       84 | T         |
|       85 | U         |
|       86 | V         |
|       87 | W         |
|       88 | X         |
|       89 | Y         |
|       90 | Z         |

## Phred Encodings

There are several different ways to encode phred scores with ascii
characters. The two most common are called phred+33 and phred+64. The
names are strange until you understand how then encoding works.

### Phred+33

To use the phred+33 encoding, take the phred quality score, add 33 to
it, then use the ascii character corresponding to the sum. For example,
using the phred+33 encoding, a quality score of 30 would be represented
with the ascii character with the ascii code of 63 (30 + 33), which is
‘?’.

### Phred+64

The phred+64 encoding works the same as the phred+33 encoding, except
you add 64 to the phred score to determine the ascii code of the quality
character. You will only find phred+64 encoding on older data, which was
sequenced several years ago. The tricky part is that there is no
indication in the FASTQ file as to which encoding was used, you have to
make an educated guess.

``` r
phred_quality %>%
  rowwise %>%
  mutate(`Phred+33`=intToUtf8(`Phred Score`+33),
         `Phred+64`=intToUtf8(`Phred Score`+64)) %>%
  kable
```

| Phred Score | Prob of Error | Phred+33 | Phred+64 |
|------------:|--------------:|:---------|:---------|
|           0 |     1.0000000 | !        | @        |
|           1 |     0.7943282 | "        | A        |
|           2 |     0.6309573 | \#       | B        |
|           3 |     0.5011872 | $        | C        |
|           4 |     0.3981072 | %        | D        |
|           5 |     0.3162278 | &        | E        |
|           6 |     0.2511886 | ’        | F        |
|           7 |     0.1995262 | (        | G        |
|           8 |     0.1584893 | )        | H        |
|           9 |     0.1258925 | \*       | I        |
|          10 |     0.1000000 | \+       | J        |
|          11 |     0.0794328 | ,        | K        |
|          12 |     0.0630957 | \-       | L        |
|          13 |     0.0501187 | .        | M        |
|          14 |     0.0398107 | /        | N        |
|          15 |     0.0316228 | 0        | O        |
|          16 |     0.0251189 | 1        | P        |
|          17 |     0.0199526 | 2        | Q        |
|          18 |     0.0158489 | 3        | R        |
|          19 |     0.0125893 | 4        | S        |
|          20 |     0.0100000 | 5        | T        |
|          21 |     0.0079433 | 6        | U        |
|          22 |     0.0063096 | 7        | V        |
|          23 |     0.0050119 | 8        | W        |
|          24 |     0.0039811 | 9        | X        |
|          25 |     0.0031623 | :        | Y        |
|          26 |     0.0025119 | ;        | Z        |
|          27 |     0.0019953 | \<       | \[       |
|          28 |     0.0015849 | =        | \\       |
|          29 |     0.0012589 | \>       | \]       |
|          30 |     0.0010000 | ?        | ^        |
|          31 |     0.0007943 | @        | \_       |
|          32 |     0.0006310 | A        | \`       |
|          33 |     0.0005012 | B        | a        |
|          34 |     0.0003981 | C        | b        |
|          35 |     0.0003162 | D        | c        |
|          36 |     0.0002512 | E        | d        |
|          37 |     0.0001995 | F        | e        |
|          38 |     0.0001585 | G        | f        |
|          39 |     0.0001259 | H        | g        |
|          40 |     0.0001000 | I        | h        |
|          41 |     0.0000794 | J        | i        |

## Why +33?

ASCII 33 is the first “normal” ASCII character that. 1 through 32
include whitespace and non-printing characters, which cannot be
identified by eye)

``` r
tibble(`ASCII #`=0:40) %>%
  rowwise %>%
  mutate(Character=intToUtf8(`ASCII #`)) %>%
  kable
```

| ASCII \# | Character |
|---------:|:----------|
|        0 |           |
|        1 |          |
|        2 |          |
|        3 |          |
|        4 |          |
|        5 |          |
|        6 |          |
|        7 |          |
|        8 |          |
|        9 |           |
|       10 |           |
|       11 |           |
|       12 |           |
|       13 |           |
|       14 |          |
|       15 |          |
|       16 |          |
|       17 |          |
|       18 |          |
|       19 |          |
|       20 |          |
|       21 |          |
|       22 |          |
|       23 |          |
|       24 |          |
|       25 |          |
|       26 |          |
|       27 |          |
|       28 |          |
|       29 |          |
|       30 |          |
|       31 |          |
|       32 |           |
|       33 | !         |
|       34 | "         |
|       35 | \#        |
|       36 | $         |
|       37 | %         |
|       38 | &         |
|       39 | ’         |
|       40 | (         |
