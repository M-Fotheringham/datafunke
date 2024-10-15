"""get_cell_counts counts cells within user-defined filters."""

from astropathdb import AstroDB
import pandas as pd
import numpy as np
from .dynamic_query import dynamic_sql
from .exprTotal import exprTotal


def get_cell_counts(
    sampleid=None,
    database=None,
    phenos=None,
    tdist_filter=None,
    rdist_filter=None,
    all_reg=False,
    excl_ln=False,
    t_hist_step=None,
    altdb=None,
    t_hist_type=None,
):

    # Preprocess inputs ========================
    if sampleid is None:
        sampleid = [114]
        print(f"Defaulting to sampleid: {sampleid}.")
    elif isinstance(sampleid, int):
        sampleid = [sampleid]

    if database is None:
        print("Defaulting to wsi02.")
        database = AstroDB(database="wsi02")

    if not isinstance(tdist_filter, tuple):
        # outer, then inner bounds
        tdist_filter = (tdist_filter, None)

    if not isinstance(rdist_filter, tuple):
        # outer, then inner bounds
        rdist_filter = (rdist_filter, None)

    if altdb is None:
        altdb = ""
    else:
        altdb = f"{altdb}."
    # ===========================================

    # Dynamic sql ==========================================================
    pheno_sql, tdist_sql, rdist_sql, ln_sql, t_hist_sql, group_sql = dynamic_sql(
        phenos, tdist_filter, rdist_filter, all_reg, excl_ln, t_hist_step, t_hist_type
    )

    # Query based on filters
    sql = f"""
    select c.sampleid, p.phenotype, c.exprphenotype,
    {t_hist_sql}count(*) c
    from {altdb}dbo.celltag c
    left join {altdb}dbo.phenotype p
        on c.ptype = p.ptype
    left join {altdb}dbo.annotations r
        on (c.sampleid = r.sampleid and r.lname = 'regression')
    left join {altdb}dbo.annotations ln
        on (c.sampleid = ln.sampleid and ln.lname = 'lymph node')
    where c.sampleid in ({",".join([f"'{x}'" for x in sampleid])})
    {tdist_sql}
    {rdist_sql}
    {ln_sql}
    {pheno_sql}
    group by c.sampleid, p.phenotype, c.exprphenotype{group_sql}
    order by c.sampleid, p.phenotype, c.exprphenotype{group_sql}   
    """
    cells = database.query(sql)
    # =======================================================================

    print(cells["tdist_bin"].max())

    # Add 0'd row for phenos with no counts
    if t_hist_step is None:
        # Group cells to identify missing exprphenotypes
        grouped = cells.groupby(["sampleid", "phenotype"])
        # List all exprs across whole cohort
        expr_list = cells["exprphenotype"].unique()

        # Where expr not present
        # Still relies on expr being present in cell df
        # Can make specific to the db
        miss_list = []
        for name, group in grouped:
            expr = group["exprphenotype"].values
            missing_bins = [e for e in expr_list if e not in expr]

            missing_rows = pd.DataFrame({"exprphenotype": missing_bins})
            missing_rows["c"] = 0
            missing_rows["sampleid"] = name[0]
            missing_rows["phenotype"] = name[1]

            miss_list.append(missing_rows)

        if miss_list:
            cells = pd.concat([cells, *miss_list])
            cells = cells.sort_values(
                ["sampleid", "phenotype", "exprphenotype"]
            ).reset_index(drop=True)

    else:
        # Group cells to identify range of tdist and missing bins
        grouped = cells.groupby(["sampleid", "phenotype", "exprphenotype"])
        # Get total range of tdist from all cells
        tdist_total = cells["tdist_bin"].unique()
        # Calculate the anticipated number of bins from overall range
        bins = np.arange(
            tdist_total.min(), tdist_total.max() + t_hist_step, t_hist_step
        )

        # Where tdist bin not present
        miss_list = []
        for name, group in grouped:
            # Get tdist bins for given group
            tdist = group["tdist_bin"].values
            # List anticipated bins not found in group tdist bins
            missing_bins = [b.astype(int) for b in bins if b not in tdist]

            # Add missing tdist bins and group labels to a df
            missing_rows = pd.DataFrame({"tdist_bin": missing_bins})
            missing_rows["c"] = 0
            missing_rows["sampleid"] = name[0]
            missing_rows["phenotype"] = name[1]
            missing_rows["exprphenotype"] = name[2]
            # Reorder columns, in case this matters (shouldn't)
            missing_rows = missing_rows[
                ["sampleid", "phenotype", "exprphenotype", "tdist_bin", "c"]
            ]

            miss_list.append(missing_rows)

        if miss_list:
            # Combine missing rows with original cells df
            # missing_df = pd.concat(miss_list, ignore_index=True)
            cells = pd.concat([cells, *miss_list]).reset_index(drop=True)

            # Convert to int (should already be)
            cells["tdist_bin"] = cells["tdist_bin"].astype(int)

    # Get a row for the total cell inclusive of all exprphenotypes
    if t_hist_step is not None:
        cells = exprTotal(cells, ["sampleid", "phenotype", "tdist_bin"])
    else:
        cells = exprTotal(cells, ["sampleid", "phenotype"])

    return cells
