import pandas as pd
import numpy as np
from astropathdb import AstroDB


def histogrammer(
    sampleid=None,
    pheno=None,
    anno=None,
    database="wsi02",
    bin_width=100,
    excl_ln=False,
    prob=False,
    samplewise=False,
):

    # bin_width in pixels

    if sampleid is None or pheno is None:
        if sampleid is None:
            sampleid = [114]
            database = "wsi02"
            print(
                f"""Sampleid not provided.
                Defaulting to {sampleid[0]} from {database}."""
            )
        elif not isinstance(
            sampleid, (list, pd.core.series.Series, np.ndarray)
        ):
            sampleid = [sampleid]
        if pheno is None:
            pheno = "CD8"
            print(f"Pheno not provided. Defaulting to {pheno}.")

    if anno is None:
        anno_col = "tdist"
    elif anno in ["tumor", "regression"]:
        anno_col = {"tumor": "tdist", "regression": "rdist"}[anno]
    else:
        print(
            """Provided anno is invalid or has no dist
            precomputed distances in celltag table. Defaulting to tumor."""
        )
        anno_col = "tdist"

    if excl_ln:
        ln_sql = "and (ln.lname is NULL or ln.ganno.STContains(ct.pos) = 0)"
    else:
        ln_sql = ""

    print(f"Querying {database}, binning {pheno} by {anno_col}.")

    db = AstroDB(database=database)

    sql = f"""
    select ct.sampleid, p.phenotype, ct.exprphenotype,
        FLOOR(ct.{anno_col} / {bin_width})*{bin_width}/2 dist_bin_um,
        COUNT(*) c
    from dbo.celltag ct
    JOIN dbo.phenotype p on ct.ptype = p.ptype
    LEFT JOIN dbo.annotations ln on ct.sampleid = ln.sampleid
        and ln.lname = 'lymph node'
    where p.phenotype = '{pheno}'
    and ct.sampleid in ({",".join(map(str, sampleid))})
    {ln_sql}
    GROUP BY ct.sampleid, p.phenotype, ct.exprphenotype,
        FLOOR(ct.{anno_col} / {bin_width})*{bin_width}/2
    """
    cells = db.query(sql)

    # Divided dist_bin by 2 to convert from pixels to um

    # Cells are assigned a dist of 32700+ in absence of annotation, exclude
    cells = cells[cells["dist_bin_um"] != 16350]

    # Get a row for the total cell inclusive of all exprphenotypes
    total = (
        cells.groupby(["sampleid", "phenotype", "dist_bin_um"], as_index=False)[
            "c"
        ]
        .sum()
        .reset_index(drop=True)
    )
    total["exprphenotype"] = "Total"

    bins = pd.concat([cells, total]).reset_index(drop=True)

    # Getting areas might be tricky, instead, can use probability
    if prob:
        if samplewise:
            bins["norm_c"] = bins["c"] / bins.groupby(
                ["sampleid", "phenotype", "exprphenotype"]
            )["c"].transform("sum")
            bins = (
                bins.groupby(["phenotype", "exprphenotype", "dist_bin_um"])[
                    "norm_c"
                ]
                .mean()
                .reset_index()
            )
        else:
            bins = (
                bins.groupby(["phenotype", "exprphenotype", "dist_bin_um"])["c"]
                .sum()
                .reset_index()
            )
            bins["prob"] = bins["c"] / bins.groupby(
                ["phenotype", "exprphenotype"]
            )["c"].transform("sum")

    return bins
