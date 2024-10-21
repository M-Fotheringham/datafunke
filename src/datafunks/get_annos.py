"""Query annotations."""

from astropathdb import AstroDB


def get_annos(
    sample=None,
    annos=None,
    database=None,
    shortcut=None,
):
    """Query annotations."""

    if sample is None:
        print("Defaulting to all samples.")
        sample_sql = ""
    else:
        if isinstance(sample, int):
            sample = [sample]

        sample_sql = f"""
        and sampleid in ({",".join([f"'{x}'" for x in sample])})
        """

    if isinstance(annos, str):
        annos = [annos]

    if database is None:
        print("Defaulting to wsi02...")
        database = AstroDB(database="wsi02")

    if shortcut is None:
        shortcut = ""

    sql = f"""
    select sampleid, a.ganno anno, a.ganno.STArea() area
    from {shortcut}dbo.annotations a
    where lname in ({",".join([f"{x}" for x in annos])})
    {sample_sql}
    """

    samples = database.query(sql)

    return samples
