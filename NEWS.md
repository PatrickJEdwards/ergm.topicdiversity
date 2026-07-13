# ergm.topicdiversity 0.2.0

- Replaced the Shannon effective-number statistic with the inverse Simpson effective number of topics: `1 / sum(p^2)`.
- Renamed output coefficients to `b1topicdiversity.inverse_simpson` and `b1topicdiversity.inverse_simpson_minus_1`.
- Updated tests to use uneven portfolios that distinguish inverse Simpson from Shannon effective diversity.
- Corrected toy network construction to use `matrix.type = "bipartite"`.
