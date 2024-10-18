"""get_areas computes tissue area from the predefined randomcell density."""

from astropathdb import AstroDB
from .dynamic_query import dynamic_sql


def get_area(
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
    """Computes tissue area from the predefined randomcell density."""

    # Preprocess inputs ========================
    if sampleid is None:
        sampleid = [114]
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
    phenos, tdist_sql, rdist_sql, ln_sql, t_hist_sql, group_sql = dynamic_sql(
        phenos,
        tdist_filter,
        rdist_filter,
        all_reg,
        excl_ln,
        t_hist_step,
        t_hist_type,
    )
    # =======================================================================

    # Normalization Unit Ratio to convert randomcell density to area
    nur = 8000.0000 / (1.004 * 1.344)

    # Query based on filters
    sql = f"""
    select c.sampleid, {t_hist_sql}count(*)/{nur} area_mm
    from {altdb}dbo.randomcell c
    left join {altdb}dbo.annotations r
        on (c.sampleid = r.sampleid and r.lname = 'regression')
    left join {altdb}dbo.annotations ln
        on (c.sampleid = ln.sampleid and ln.lname = 'lymph node')
    where c.sampleid in ({",".join([f"'{x}'" for x in sampleid])})
    {tdist_sql}
    {rdist_sql}
    {ln_sql}
    group by c.sampleid{group_sql}
    """
    area = database.query(sql)

    if t_hist_step is not None:
        area = area.sort_values(["sampleid", "tdist_bin"]).reset_index(
            drop=True
        )
    else:
        area = area.sort_values("sampleid").reset_index(drop=True)

    return area
