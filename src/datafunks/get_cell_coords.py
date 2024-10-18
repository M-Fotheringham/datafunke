"""Get cell coordinates for a specified phenotype and analysis boundary"""

from astropathdb import AstroDB


def get_cell_coords(
    sampleid=None,
    pheno=None,
    database="wsi02",
    tdist_filter=None,
    excl_ln=False,
    reg_only=False,
    all_reg=False,
):
    """Get cell coordinates for a specified phenotype and analysis boundary"""

    if sampleid is None or pheno is None:
        if sampleid is None:
            sampleid = 114
            database = "wsi02"
            print(
                f"Sampleid not provided. Defaulting to {114} from {database}."
            )
        if pheno is None:
            pheno = "CD8"
            print(f"Pheno not provided. Defaulting to {pheno}.")

    print(f"Querying {database}, {sampleid} {pheno}.")

    db = AstroDB(database=database)

    if not isinstance(tdist_filter, tuple):
        # outer, then inner bounds
        tdist_filter = (tdist_filter, None)

    if tdist_filter[0] is not None and tdist_filter[1] is not None:
        if all_reg:
            tdist_filt = f"""
            and ((ct.tdist <= {tdist_filter[0]}
            and ct.tdist >= {tdist_filter[1]})
                or d.ganno.STContains(ct.pos) = 1)
            """
        else:
            tdist_filt = f"""
            and (ct.tdist <= {tdist_filter[0]}
            and ct.tdist >= {tdist_filter[1]})
            """
    elif tdist_filter[0] is not None:
        if all_reg:
            tdist_filt = f"""
            and (ct.tdist <= {tdist_filter[0]}
                or d.ganno.STContains(ct.pos) = 1)
            """
        else:
            tdist_filt = f"""
            and ct.tdist <= {tdist_filter[0]}
            """
    elif tdist_filter[1] is not None:
        if all_reg:
            tdist_filt = f"""
            and (ct.tdist >= {tdist_filter[1]}
                or d.ganno.STContains(ct.pos) = 1)
            """
        else:
            tdist_filt = f"""
            and ct.tdist >= {tdist_filter[1]}
            """
    else:
        tdist_filt = ""

    if excl_ln:
        # Exclude cells in lymph node for samples with ln annotations
        ln_sql = """
        and (c.ganno.STContains(ct.pos) = 0 or c.lname is NULL)
        """
    else:
        ln_sql = ""

    if reg_only:
        # Only cells in regression
        reg_sql = """
        and (d.ganno.STContains(ct.pos) = 1 and d.lname = 'regression')
        """
    else:
        reg_sql = ""

    sql = f"""
    select p.phenotype, ct.exprphenotype, ct.tdist, ct.px, ct.py
    from dbo.celltag ct
    JOIN dbo.phenotype p on ct.ptype = p.ptype
    LEFT JOIN dbo.annotations c on ct.sampleid = c.sampleid
        and c.lname = 'lymph node'
    LEFT JOIN dbo.annotations d on ct.sampleid = d.sampleid
        and d.lname = 'regression'
    where p.phenotype = '{pheno}'
    and ct.sampleid = {sampleid}
    {ln_sql}
    {reg_sql}
    {tdist_filt}
    """
    cells = db.query(sql)

    return cells
