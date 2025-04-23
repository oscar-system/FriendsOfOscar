# Wronski Pairs of Honeycomb Curves

This is code accompanying the [paper of the same
name](https://arxiv.org/abs/2411.10776) by [Laura
Casabella](https://sites.google.com/view/lauracasabella), [Michael
Joswig](https://page.math.tu-berlin.de/~joswig/) and [Rafael
Mohr](https://mathexp.eu/mohr/).

# Usage

Run julia with the project defined in this directory.
That is, after `cd` into this directory, run the following

```bash
julia  --project
```

After then calling `using Honeycomb`, the data in Table 1 of the
abovementioned paper can be reproduced by running the function
`record_data_numeric(17, "<l_func>")` where `<l_func>` is `rho` for
the left column and `min` for the second column.

The data in Table 2 can be reproduced by running the function
`record_data_solve(9, "<l_func>")` in a similar manner and the data
in Table 3 can be reproduced using `record_data(11, "<l_func>")`.
