"""Get  the cell coordinates for a specified phenotype and analysis boundary"""

import pandas as pd
import geopandas as gpd
from shapely.wkt import loads
from astropathdb import AstroDB
from .anno_check import anno_check


def get_cell_den(
    sampleid=None,
    pheno=None,
    database="wsi02",
    tdist_filter=None,
    excl_ln=False,
    reg_only=False,
    all_reg=False,
):

    print(f"Querying {database}...")

    if database == "wsi34":
        database = "wsi02"
        wsi34_shortcut = "wsi34."
        print("Using wsi02.wsi34 to circumvent access issue...")
    else:
        wsi34_shortcut = ""

    db = AstroDB(database=database)

    if sampleid is None:
        if reg_only:
            sampleid_sql = f"""
            ct.sampleid in ({",".join(map(str, anno_check("regression", db, shortcut=wsi34_shortcut)))})
            """
            print(
                f"Sampleid not provided. Defaulting to all samples with regression from {database}."
            )
        else:
            sampleid_sql = f"""
            ct.sampleid in ({",".join(map(str, anno_check("good tissue", db, shortcut=wsi34_shortcut)))})
            """
            print(f"Sampleid not provided. Defaulting to all samples from {database}.")
    else:
        if not isinstance(sampleid, list):
            sampleid = [sampleid]
        sampleid_sql = f"""
        ct.sampleid in ({",".join(map(str, sampleid))})
        """

    if pheno is None:
        pheno_sql = ""
        print(f"Pheno not provided. Defaulting to all phenotypes.")
    else:
        if not isinstance(pheno, list):
            pheno = [pheno]
        pheno_sql = f"""
        and p.phenotype in ({",".join([f"'{x}'" for x in pheno])})
        """

    if not isinstance(tdist_filter, tuple):
        # outer, then inner bounds
        tdist_filter = (tdist_filter, None)

    if tdist_filter[0] is not None and tdist_filter[1] is not None:
        if all_reg:
            tdist_sql = f"""
            and ((ct.tdist <= {tdist_filter[0]}
            and ct.tdist >= {tdist_filter[1]}) or (d.lname is not NULL and d.ganno.STContains(ct.pos) = 1))
            """
        else:
            tdist_sql = f"""
            and (ct.tdist <= {tdist_filter[0]}
            and ct.tdist >= {tdist_filter[1]})
            """
        t_anno_sql = f"""
        a.ganno.STBuffer({tdist_filter[0]}).STDifference(a.ganno.STBuffer({tdist_filter[1]})).STIntersection(b.ganno)
        """
    elif tdist_filter[0] is not None:
        if all_reg:
            tdist_sql = f"""
            and (ct.tdist <= {tdist_filter[0]} or (d.lname is not NULL and d.ganno.STContains(ct.pos) = 1))
            """
        else:
            tdist_sql = f"""
            and ct.tdist <= {tdist_filter[0]}
            """
        t_anno_sql = f"""
        a.ganno.STBuffer({tdist_filter[0]}).STIntersection(b.ganno)
        """
    elif tdist_filter[1] is not None:
        if all_reg:
            tdist_sql = f"""
            and (ct.tdist >= {tdist_filter[1]} or (d.lname is not NULL and d.ganno.STContains(ct.pos) = 1))
            """
        else:
            tdist_sql = f"""
            and ct.tdist >= {tdist_filter[1]}
            """
        t_anno_sql = f"""
        b.ganno.STDifference(a.ganno.STBuffer({tdist_filter[1]}))
        """
    else:
        tdist_sql = ""
        t_anno_sql = "b.ganno"

    if excl_ln:
        ln_sql = f"""
        and (ln.ganno.STContains(ct.pos) = 0 or ln.lname is NULL)
        """
    else:
        ln_sql = ""

    if reg_only:
        reg_sql = f"""
        and d.lname = 'regression'
        and d.ganno.STContains(ct.pos) = 1
        """
    else:
        reg_sql = ""

    print("Counting cells...")

    sql = f"""
    SELECT ct.sampleid, p.phenotype, ct.exprphenotype, COUNT(*) c
    FROM {wsi34_shortcut}dbo.celltag ct
    JOIN {wsi34_shortcut}dbo.phenotype p on ct.ptype = p.ptype
    LEFT JOIN {wsi34_shortcut}dbo.annotations d on ct.sampleid = d.sampleid and d.lname = 'regression'
    LEFT JOIN {wsi34_shortcut}dbo.annotations ln on ct.sampleid = ln.sampleid and ln.lname = 'lymph node'
    WHERE {sampleid_sql}
    {reg_sql}
    {ln_sql}
    {pheno_sql}
    {tdist_sql}
    GROUP BY ct.sampleid, p.phenotype, ct.exprphenotype
    """
    cells = db.query(sql)

    if database == "wsi02":
        cells = cells[~((cells["sampleid"].isin([547, 566])) & (cells["phenotype"] == "FoxP3CD8"))]

    samples = cells["sampleid"].unique()

    # Add null value to avoid syntax error in annotation sql query
    if len(samples) == 0:
        samples = ["NULL"]

    if excl_ln:
        # Exclude cells in lymph node for samples with ln annotations
        ln_anno_sql = ".STDifference(c.ganno)"
    else:
        ln_anno_sql = ""

    if reg_only:
        # Include only cells in regression for samples with reg annotations
        reg_anno_sql = ".STIntersection(d.ganno)"
    else:
        reg_anno_sql = ""

    if all_reg:
        # Include all regression regardless of tdist filter
        all_reg_sql = ".STUnion(d.ganno).STIntersection(b.ganno)"
        all_reg_sql2 = ".STIntersection(d.ganno)"
    else:
        all_reg_sql = ""
        all_reg_sql2 = ""

    print("Calculating areas...")

    areas_list = []
    for sample in samples:

        print(sample)

        sql = f"""
        select b.sampleid,
        CASE
            when a.lname is not NULL and d.lname is not NULL and c.lname is not NULL then {t_anno_sql}{reg_anno_sql}{all_reg_sql}{ln_anno_sql}
            when a.lname is not NULL and d.lname is not NULL and c.lname is NULL then {t_anno_sql}{reg_anno_sql}{all_reg_sql}
            when a.lname is not NULL and d.lname is NULL and c.lname is not NULL then {t_anno_sql}{ln_anno_sql}
            when a.lname is not NULL and d.lname is NULL and c.lname is NULL then {t_anno_sql}
            when a.lname is NULL and c.lname is not NULL then b.ganno{reg_anno_sql}{all_reg_sql2}{ln_anno_sql}
            when a.lname is NULL and c.lname is NULL then b.ganno{reg_anno_sql}{all_reg_sql2}
        END anno
        from {wsi34_shortcut}dbo.annotations b
        left join {wsi34_shortcut}dbo.annotations a on b.sampleid = a.sampleid and a.lname = 'tumor'
        left join {wsi34_shortcut}dbo.annotations c on b.sampleid = c.sampleid and c.lname = 'lymph node'
        left join {wsi34_shortcut}dbo.annotations d on b.sampleid = d.sampleid and d.lname = 'regression'
        where b.sampleid in ({sample})
        and b.lname = 'good tissue'
        """
        # ",".join(map(str, samples))
        sample_area = db.query(sql)
        areas_list.append(sample_area)

    areas = pd.concat(areas_list)

    # Add something here to fill in empty (dropped) rows

    # It's faster to calculate areas outside of the query for some reason
    areas["area"] = [gpd.GeoSeries(loads(x)).area[0] * 2.5e-7 for x in areas["anno"]]

    # Combine with cells to calculate density
    den = pd.merge(cells, areas, on="sampleid").reset_index(drop=True)

    # In case of providing plots
    # den["pos"] = den["pos"].apply(loads).apply(gpd.GeoSeries)

    # Get a row for the total cell inclusive of all exprphenotypes
    total = (
        den.groupby(["sampleid", "phenotype", "anno", "area"], as_index=False)["c"]
        .sum()
        .reset_index(drop=True)
    )
    total["exprphenotype"] = "Total"

    den = pd.concat([den, total]).reset_index(drop=True)

    den["density"] = den["c"] / den["area"]

    den["filters"] = (
        f"tdist:{tdist_filter}_all_reg:{all_reg}_reg_only:{reg_only}_excl_ln:{excl_ln}"
    )

    print("Done.")

    return den
