"""Get Total row for queries with exprPhenotype."""

import pandas as pd


def expr_total(df, cols, sum_col="c", expr_col="exprphenotype"):
    """Get Total row for queries with exprPhenotype.
    df: DataFrame
    cols: list of columns containing grouping variables
    sum_col: column containing cell counts to group
    expr_col: column in which to add Total label
    """

    # Create a 'total' cell row inclusive of all exprphenotypes per lineage
    total = df.groupby(cols, as_index=False)
    total = total[sum_col].sum().reset_index(drop=True)

    total[expr_col] = "Total"  # label total

    # Combine with original df and sort
    df_total = pd.concat([df, total]).reset_index(drop=True)

    df_total = df_total.sort_values(cols).reset_index(drop=True)

    return df_total
