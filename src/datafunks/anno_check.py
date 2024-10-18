# Check which cases is in the current db have certain annotations
def anno_check(annos_to_check, database, shortcut=None):

    if isinstance(annos_to_check, str):
        annos_to_check = [annos_to_check]

    if shortcut is None:
        shortcut = ""

    sql = f"""
    select distinct sampleid
    from {shortcut}dbo.annotations
    where lname in ({",".join([f"'{x}'" for x in annos_to_check])})
    """

    samples = database.query(sql)["sampleid"].to_list()

    return samples
