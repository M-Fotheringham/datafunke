"""Get Total row for queries with exprPhenotype."""

import pandas as pd


def expr_total(df, cols, sum_col="c", expr_col="exprphenotype"):
    """Get Total row for queries with exprPhenotype."""

    # Get a row for the total cell inclusive of all exprphenotypes
    total = df.groupby(cols, as_index=False)
    total = total[sum_col].sum().reset_index(drop=True)

    total[expr_col] = "Total"

    df_total = pd.concat([df, total]).reset_index(drop=True)

    df_total = df_total.sort_values(cols).reset_index(drop=True)

    # Add 0 for skipped (empty) rows at some point

    return df_total
