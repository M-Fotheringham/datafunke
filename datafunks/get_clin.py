"""Queries the clinical data for the relevant cohort and matches by patient.
Of little use- all clin tables are different and wsi02 is not up-to-date."""

from astropathdb import AstroDB


def get_clin(sampleid=None, database="wsi02"):

    if database == "wsi02":
        print("Recall that wsi02 clin info not up-to-date in database.")

    if sampleid is not None:
        if not isinstance(sampleid, list):
            sampleid = [sampleid]
        sampleid_sql = f"""where sampleid in ({",".join(map(str, sampleid))})"""
    else:
        print(f"Sampleid not provided, defaulting to all in {database}.")
        sampleid_sql = ""

    if database in ["wsi06", "wsi02"] and sampleid is not None:
        sampleid_sql = f"""
        where slideid in (
            select
                CASE
                    WHEN CHARINDEX('_', slideid) > 0
                    THEN SUBSTRING(slideid, 1, CHARINDEX('_', slideid) - 1)
                    ELSE slideid
                END slideid
            from dbo.samples
            where sampleid in ({",".join(map(str, sampleid))})
            )
        """

    db = AstroDB(database=database)

    sql = f"""
    select *
    from dbo.clinical c
    {sampleid_sql}
    """
    clin = db.query(sql)

    return clin
