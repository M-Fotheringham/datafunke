def dynamic_sql(
    phenos=None,
    tdist_filter=None,
    rdist_filter=None,
    all_reg=False,
    excl_ln=False,
    t_hist_step=None,
    t_hist_type=None,
):
    """This saves some tedious repition in sql queries for common tasks."""

    if (t_hist_step is not None) & (t_hist_type != "fractional reg"):
        # convert to pixels
        t_hist_step = t_hist_step * 2

    # Phenotype filter
    if phenos is not None:
        if isinstance(phenos, str):
            phenos = [phenos]
        pheno_sql = f"""
        and p.phenotype in ({",".join([f"'{x}'" for x in phenos])})
        """
    else:
        pheno_sql = ""

    if all_reg:
        # Dynamic tdist filtering with all regression
        if (tdist_filter[0] is not None) and (tdist_filter[1] is not None):
            tdist_sql = f"""
            and ((tdist <= {tdist_filter[0]} and tdist > {tdist_filter[1]})
                or rdist <= 0)
            """
        elif tdist_filter[0] is not None:
            tdist_sql = f"""
            and (tdist <= {tdist_filter[0]} or rdist <= 0)
            """
        elif tdist_filter[0] is not None:
            tdist_sql = f"""
            and (tdist > {tdist_filter[1]} or rdist <= 0)
            """
        else:
            tdist_sql = ""

    else:
        # Dynamic tdist filtering without all regression
        if (tdist_filter[0] is not None) and (tdist_filter[1] is not None):
            tdist_sql = f"""
            and (tdist <= {tdist_filter[0]} and tdist > {tdist_filter[1]})
            """
        elif tdist_filter[0] is not None:
            tdist_sql = f"""
            and tdist <= {tdist_filter[0]}
            """
        elif tdist_filter[0] is not None:
            tdist_sql = f"""
            and tdist > {tdist_filter[1]}
            """
        else:
            tdist_sql = ""

    # Dynamic rdist filtering
    if (rdist_filter[0] is not None) and (rdist_filter[1] is not None):
        rdist_sql = f"""
        and (rdist <= {rdist_filter[0]} and rdist > {rdist_filter[1]})
        """
    elif rdist_filter[0] is not None:
        rdist_sql = f"""
        and rdist <= {rdist_filter[0]}
        """
    elif rdist_filter[0] is not None:
        rdist_sql = f"""
        and rdist > {rdist_filter[1]}
        """
    else:
        rdist_sql = ""

    # Exclude lymph node
    if excl_ln:
        ln_sql = "and (ln.lname is NULL or ln.ganno.STContains(c.pos) = 0)"
    else:
        ln_sql = ""

    # Bin tdist to create histogram of areas, convert tdist to um
    if t_hist_type == "fractional reg":
        t_hist_sql = f"""
        floor((tdist/(tdist - rdist))*100/
        {t_hist_step})*{t_hist_step} tdist_bin,
        """
        group_sql = f"""
                     , floor((tdist/(tdist - rdist))*100/
                     {t_hist_step})*{t_hist_step}
                    """
    elif t_hist_step is not None:
        t_hist_sql = f"""
        floor(tdist/{t_hist_step})*{t_hist_step}/2 tdist_bin,
        """
        group_sql = f""", floor(tdist/{t_hist_step})*{t_hist_step}/2"""
    else:
        t_hist_sql = ""
        group_sql = ""

    return pheno_sql, tdist_sql, rdist_sql, ln_sql, t_hist_sql, group_sql
