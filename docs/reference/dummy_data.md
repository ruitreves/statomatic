# Produces a dummy data set from a normal distribution and sample_info object.

Produces a dummy data set from a normal distribution and sample_info
object.

## Usage

``` r
dummy_data(
  group_min = 2,
  group_max = 20,
  dim_min = 1,
  dim_max = 50,
  param_min = 0,
  param_max = 20
)
```

## Arguments

- group_min:

  Minimum random number of groups. Defaults to 3

- group_max:

  Maximum random number of groups. Defaults to 20. Note: group_max \*
  dim max must be \<= 1000

- dim_min:

  Minimum random dimension of data frame. Defaults to 1

- dim_max:

  Maximum random dimension of data frame. Defaults to 50. Note:
  group_max \* dim max must be \<= 1000

- param_min:

  Minimum random parameters for normal distribution (i.e., mean and
  variance). Defaults to 0

- param_max:

  Maximum random parameters for normal distribution. Defaults to 20

## Value

A list, element one is the randomly generated data frame, element two is
sample_info. Access with double brackets: `[[1]]`, `[[2]]`.
