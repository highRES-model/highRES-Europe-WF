import calendar
import datetime
import itertools
import re
import warnings

import numpy as np
import pandas as pd
import pathlib


def wrapdd(d, out_par, otype, outfile=""):
    d = np.atleast_2d(d)

    if otype == "set":
        top = np.array([["set"], [out_par + " /"]])
        bottom = np.array([["/"], [""]])

    if otype == "scalar":
        top = np.array([["scalar"], [out_par + " /"]])
        bottom = np.array([["/"], [""]])

    if otype == "parameter":
        top = np.array([["parameter", ""], [out_par + " /", ""]])
        bottom = np.array([["/", ""], ["", ""]])

    d = np.concatenate((top, d, bottom), axis=0)

    if outfile != "":
        np.savetxt(outfile, d, delimiter=" ", fmt="%s")
    else:
        return d


def data2dd(data, sets, outfile="", all_combin=False, rounddp=8):
    # rounddp=5

    # each set in sets needs to be 1D at the moment

    # sets and data must be in correct order -> last set must be column
    # headers, previous sets are rows

    if type(data) != "numpy.ndarray":
        data = np.array(data)
    sets = np.array(sets, dtype="object")

    if len(sets) > 0:
        data = np.round(data.astype(float), rounddp)
        if len(sets) == 1:
            sets_out = sets[0]
            sets_new = np.array(sets_out).astype(str)

            data_out = np.hstack((sets_new.reshape(-1, 1), data.reshape(-1, 1)))
        else:
            if all_combin:
                sets_out = [item for item in list(itertools.product(*sets))]

                sets_out = np.char.array(sets_out)

                sets_new = sets_out[:, 0]

                for i in range(sets_out.shape[1] - 1):
                    sets_new = sets_new + "." + sets_out[:, i + 1]

                data_out = np.hstack((sets_new.reshape(-1, 1), data.reshape(-1, 1)))

            else:
                lens = np.array([item.shape[0] for item in sets])

                sets_new = []
                for i, s in enumerate(sets):
                    if lens[i] != lens.max():
                        sets_new.append(np.repeat(s, lens.max()))
                    else:
                        sets_new.append(s)

                sets_new = np.array(
                    [".".join(item) for item in np.char.array(list(zip(*sets_new)))]
                )

                data_out = np.hstack((sets_new.reshape(-1, 1), data.reshape(-1, 1)))

    else:
        data_out = data.reshape(-1, 1)
    if outfile != "":
        np.savetxt(outfile, data_out, fmt="%s", delimiter=" ")
    else:
        return data_out


def temporal2dd(dstart, dend, opath, temporaloutputpath):
    # nyears=(dend.year-dstart.year)+1
    years = np.arange(dstart.year, dend.year + 1)
    ntime = ((dend - dstart).total_seconds() / 3600) + 1

    hr2yr = []
    for nyr, yr in enumerate(years):
        shour = ((datetime.datetime(yr, 1, 1, 0) - dstart).total_seconds() / 3600) + 1
        if dend.year == yr:
            ehour = ((dend - dstart).total_seconds() / 3600) + 1
        else:
            ehour = (
                (datetime.datetime(yr, 12, 31, 23) - dstart).total_seconds() / 3600
            ) + 1

        hrs = np.arange(shour - 1, ehour).astype(int)
        hr2yr.append(list(zip(np.repeat(nyr, hrs.shape[0]).astype(int), hrs)))

    hr2yr = np.char.array(np.vstack(hr2yr).astype(str))

    hr2yr = (hr2yr[:, 0] + "." + hr2yr[:, 1]).reshape(-1, 1)

    if dstart.year != dend.year:
        yr = str(dstart.year) + "-" + str(dend.year)
    else:
        yr = str(dstart.year)

    out = []

    out.append(wrapdd(np.arange(ntime).reshape(-1, 1).astype(int), "h", "set"))
    out.append(
        wrapdd(np.arange(years.shape[0]).reshape(-1, 1).astype(int), "yr", "set")
    )
    out.append(wrapdd(hr2yr, "hr2yr_map", "set"))

    np.savetxt(temporaloutputpath, np.concatenate(out, axis=0), fmt="%s")


def trans_links(root, f, aggregated_regions, out="work"):
    tech_type = "trans"

    links_allowed = pd.read_excel(
        f, sheet_name="transmission_allowed", skiprows=1, engine="calamine"
    ).query("Zone1 == @aggregated_regions and Zone2 == @aggregated_regions")

    links_out = np.array(
        links_allowed["Zone1"]
        + "."
        + links_allowed["Zone2"]
        + "."
        + links_allowed["Tech"]
    )

    links_tech = pd.unique(links_allowed["Tech"])

    set_outdd = []
    param_outdd = []

    set_outdd.append(wrapdd(data2dd(links_tech, []), "trans", "set"))
    set_outdd.append(wrapdd(data2dd(links_out, []), "trans_links", "set"))

    # data2dd(links_out,[],outfile=root/out/(tech_type+"_links.dd"))

    out_data = ["links_dist", "links_cap","links_lim_cap"]

    for od in out_data:
        out_par = tech_type + "_" + od

        param_outdd.append(
            wrapdd(data2dd(links_allowed[od].values, [links_out]), out_par, "parameter")
        )

    params = pd.read_excel(f, sheet_name="transmission", skiprows=1, engine="calamine")

    nans = params.isnull()
    params = params.where(~nans, other=0.0)

    out_data = ["loss", "line_capex", "sub_capex", "varom"]

    #    trans_set=params["Technology Name (highRES)"]

    params = params[params["Technology Name (highRES)"].isin(links_tech)]

    for od in out_data:
        out_par = tech_type + "_" + od

        param_outdd.append(
            wrapdd(
                data2dd(params[od], [params["Technology Name (highRES)"]]),
                out_par,
                "parameter",
            )
        )

    param_outdd = np.concatenate(param_outdd, axis=0)
    set_outdd = np.concatenate(set_outdd, axis=0)

    pad = np.repeat(np.array(""), set_outdd.shape[0]).reshape(set_outdd.shape[0], 1)

    outdd = np.concatenate((np.hstack((set_outdd, pad)), param_outdd), axis=0)

    #    outdd=np.concatenate(outdd,axis=0)

    np.savetxt(root / out / (tech_type + ".dd"), outdd, delimiter=" ", fmt="%s")
    
def add_vre_connection_costs(
            root,
            out,
            f_techno,
            psys_scen,
            vre_connection_dists
            ):
    
    dists=pd.read_csv(vre_connection_dists,index_col=0)
    
    capex=pd.read_excel(f_techno, 
                             sheet_name="transmission", 
                             skiprows=1, 
                             engine="calamine")
    
    conn_costs=(dists.div(100)
            .mul(capex.loc[capex["Technology Name (highRES)"]=="HVDC_Windoffshore","line_capex"].values[0]))
    
    conn_costs.loc[(dists["dist"]>0.0),"dist"]+=capex.loc[capex["Technology Name (highRES)"]=="HVDC_Windoffshore","sub_capex"].values[0]
    
    conn_costs=pd.DataFrame(wrapdd(conn_costs.round(1).reset_index().values,
                 "gen_vreconnection","parameter"))
    
    conn_costs["out"]=conn_costs[0]+" "+conn_costs[1].astype(str)
    
    d=pd.read_csv(root / out / (psys_scen+"_gen.dd"),skip_blank_lines=False,header=None)
    
    d_out=np.concatenate(
        (d.values,
        conn_costs.loc[:,["out"]].values),
        axis=0)
    
    np.savetxt(
            root / out / (psys_scen+"_gen.dd"), d_out, delimiter=" ", fmt="%s"
        )

def co2target2dd(co2targets_db, 
                 co2target_out, 
                 esys_scen, 
                 co2_target_type,
                 co2_target_extent,
                 planning_horizon):
    
    d = pd.read_csv(co2targets_db)
    d["year"] = d["year"].astype(str).str.strip()
    dout = (d
            .query("(case == @esys_scen) \
                   and (type == @co2_target_type) \
                   and (extent == @co2_target_extent) \
                   and (year == @planning_horizon)"))

                   
    if co2_target_extent == "all":
        wrapdd(dout["target"].values[0], "co2_target", "scalar", outfile=co2target_out)
    if co2_target_extent == "zonal":
        wrapdd(data2dd(dout["target"].values,
                       [dout["zone"].values])
               ,"co2_target", "parameter", outfile=co2target_out)
        

def getzlims(lim, techs, zones):
    lim = lim.loc[(lim["Year"] == 2030) & (lim["Technology"].isin(techs)), :]

    if lim.empty:
        return np.array([])

    have_lim = ~lim.loc[:, "limtype"].isnull()

    if (~have_lim).any():
        warnings.warn(
            "Warning, missing zonal new capacity limits for:"
            + ", ".join(lim.loc[~have_lim, "Technology"])
        )

    lim = lim.loc[have_lim, :]
    para_lim = lim["parameter"].drop_duplicates()

    outlims = []

    for p_lim in para_lim:
        nl = []
        for _, row in lim.loc[lim["parameter"] == p_lim, :].iterrows():
            zones = row[zones].index.values
            limval = row[zones].values.astype(float) / 1e3
            limval[np.isnan(limval)] = 0.0

            tech = np.atleast_1d(row["Technology"])
            limtype = np.atleast_1d(row["limtype"])

            nl.append(data2dd(limval.T, [zones, tech, limtype]))

        outlims.append(wrapdd(np.concatenate(nl, axis=0), p_lim, "parameter"))

    return np.concatenate(outlims, axis=0)


def getrlims(lim, techs, zones, exist_agg):
    # TODO capcity units are fixed to GW here, need to add flexibility

    lim = lim.loc[
        (lim["Year"] == 2030)
        & (lim["Technology"].isin(techs))
        & (lim["zone"].isin(zones)),
        :,
    ]

    if lim.empty:
        return np.array([])

    have_lim = ~lim.loc[:, "limtype"].isnull()
    lim = lim.loc[have_lim, :]
    para_lim = lim["parameter"].drop_duplicates()

    outlims = []

    for p_lim in para_lim:
        limval = lim.loc[(lim["parameter"] == p_lim), :]

        if exist_agg == "nuts2":
            nl = data2dd(
                limval["value"] / 1e3,
                [
                    limval["Technology"],
                    limval["zone"],
                    limval["region"],
                    limval["limtype"],
                ],
            )

        if exist_agg == "region":
            #limval = limval.groupby(["zone"], as_index=False).agg(
            # to consider the other technologies
            limval = limval.groupby(["Technology","zone"], as_index=False).agg(
                {
                    #"Technology": "first", # removed because is already in groupby
                    "Year": "first",
                    "parameter": "first",
                    "limtype": "first",
                    "value": "sum",
                }
            )

            nl = data2dd(
                limval["value"] / 1e3,
                [
                    limval["Technology"],
                    limval["zone"],
                    limval["zone"],
                    limval["limtype"],
                ],
            )

        outlims.append(wrapdd(nl, p_lim, "parameter"))

    return np.concatenate(outlims, axis=0)


def read_dd_blocks(f):
    lines = pathlib.Path(f).read_text().splitlines()
    blocks, order = {}, []
    i, n = 0, len(lines)
    while i < n:
        line = lines[i].strip()
        if line in ("parameter", "set", "scalar"):
            block_type = line
            i += 1
            name = lines[i].strip().split("/")[0].strip()
            i += 1
            rows = []
            while i < n and lines[i].strip() != "/":
                parts = lines[i].split()
                if len(parts) == 2:
                    rows.append((parts[0], parts[1]))
                elif len(parts) == 1:
                    rows.append((parts[0], None))
                i += 1
            i += 1
            blocks[name] = {"type": block_type, "rows": rows}
            order.append(name)
        else:
            i += 1
    return blocks, order

def write_dd_blocks(blocks, order, outfile):
    lines = []
    for name in order:
        block = blocks[name]
        lines.append(block["type"])
        lines.append(name + " /")
        for key, value in block["rows"]:
            lines.append(key if value is None else f"{key} {value}")
        lines.append("/")
    pathlib.Path(outfile).write_text("\n".join(lines) + "\n")


def patch_exist_block(blocks, block_name, ledger_slice, key_cols):
    original_rows = blocks[block_name]["rows"]
    tech_pos = key_cols.index("Technology")
    limtype_lookup = {}
    for key, _ in original_rows:
        parts = key.split(".")
        limtype_lookup[parts[tech_pos]] = parts[-1]

    if ledger_slice.empty:
        blocks[block_name]["rows"] = []
        return

    agg = ledger_slice.groupby(key_cols, as_index=False)["capacity_mw"].sum()
    new_rows = []
    for _, row in agg.iterrows():
        key_parts = [str(row[c]) for c in key_cols]
        limtype = limtype_lookup.get(row["Technology"], "FX")
        key = ".".join(key_parts) + "." + limtype
        value = row["capacity_mw"] / 1000
        new_rows.append((key, str(value)))
    blocks[block_name]["rows"] = new_rows

def apply_capacity_retirement(gen_ddfile, store_ddfile, ledger_path, f_techno, planning_horizon, spatials):
    ledger = pd.read_csv(ledger_path)
    survivors = ledger.loc[ledger["installed_year"] + ledger["lifetime"] > planning_horizon, :]

    gen_techs = set(pd.read_excel(f_techno, sheet_name="gen", skiprows=1, engine="calamine")["Technology Name (highRES)"])
    store_techs = set(pd.read_excel(f_techno, sheet_name="store", skiprows=1, engine="calamine")["Technology Name (highRES)"])

    gen_pcap = survivors.loc[(survivors["capacity_type"] == "pcap") & survivors["Technology"].isin(gen_techs), :]
    gen_pcap_z_direct = gen_pcap.loc[gen_pcap["region"].isna(), :]
    gen_pcap_r_native = gen_pcap.loc[gen_pcap["region"].notna(), :]

    gen_pcap_z_from_r = gen_pcap_r_native.groupby(["zone", "Technology"], as_index=False)["capacity_mw"].sum()
    gen_pcap_z = pd.concat([gen_pcap_z_direct, gen_pcap_z_from_r], ignore_index=True)

    if spatials == "region":
        gen_pcap_r = gen_pcap_r_native.groupby(["Technology", "zone"], as_index=False)["capacity_mw"].sum()
        gen_pcap_r["region"] = gen_pcap_r["zone"]
    else:
        gen_pcap_r = gen_pcap_r_native

    gen_ecap_z = survivors.loc[
        (survivors["capacity_type"] == "ecap") & survivors["Technology"].isin(gen_techs) & survivors["region"].isna(), :
    ]

    store_pcap_z = survivors.loc[(survivors["capacity_type"] == "pcap") & survivors["Technology"].isin(store_techs), :]
    store_ecap_z = survivors.loc[(survivors["capacity_type"] == "ecap") & survivors["Technology"].isin(store_techs), :]

    gen_blocks, gen_order = read_dd_blocks(gen_ddfile)
    patch_exist_block(gen_blocks, "gen_exist_pcap_z", gen_pcap_z, ["zone", "Technology"])
    patch_exist_block(gen_blocks, "gen_exist_ecap_z", gen_ecap_z, ["zone", "Technology"])
    patch_exist_block(gen_blocks, "gen_exist_pcap_r", gen_pcap_r, ["Technology", "zone", "region"])
    write_dd_blocks(gen_blocks, gen_order, gen_ddfile)

    store_blocks, store_order = read_dd_blocks(store_ddfile)
    patch_exist_block(store_blocks, "store_exist_pcap_z", store_pcap_z, ["zone", "Technology"])
    patch_exist_block(store_blocks, "store_exist_ecap_z", store_ecap_z, ["zone", "Technology"])
    write_dd_blocks(store_blocks, store_order, store_ddfile)

def scen2dd(
    co2budgetddlocation,
    root,
    fin,
    params2write,
    scen_db,
    run,
    esys,
    zones,
    out="work",
    esys_cap=False,
    exist_cap=False,
    exist_agg="region",
):
    # co2lim2dd(co2budgetddlocation, root, run, esys, scen_db, out=out)

    scen = pd.read_excel(scen_db, sheet_name="scenario_tech_definition", skiprows=0)

    scen = scen.loc[(scen["Psys Scenario"] == run), :]

    techs = scen["Technology Name (highRES)"]
    tech_class = ["gen", "store"]

    for tech_type in tech_class:
        set_outdd = []
        param_outdd = []

        params = pd.read_excel(fin, sheet_name=None, skiprows=1, engine="calamine")[
            tech_type
        ]

        params = params[params["Technology Name (highRES)"].isin(techs)]

        sets = params2write[tech_type]["set"]
        data_out = params2write[tech_type]["parameter"]

        params = params.merge(scen, on="Technology Name (highRES)")

        # new_lim=pd.read_excel(fin,sheet_name=tech_type+"_lim_z",skiprows=0)

        param_outdd.append(
            getzlims(
                pd.read_excel(
                    fin, sheet_name=tech_type + "_lim_z", skiprows=0, engine="calamine"
                ),
                techs,
                zones,
            )
        )

        if exist_cap:
            zlims = getzlims(
                pd.read_excel(
                    fin,
                    sheet_name=tech_type + "_exist_z",
                    skiprows=0,
                    engine="calamine",
                ),
                techs,
                zones,
            )

            if zlims.size != 0:
                param_outdd.append(zlims)

            if tech_type == "gen":
                rlims = getrlims(
                    pd.read_excel(
                        fin,
                        sheet_name=tech_type + "_exist_r",
                        skiprows=0,
                        engine="calamine",
                    ),
                    techs,
                    zones,
                    exist_agg,
                )

                if rlims.size != 0:
                    param_outdd.append(rlims)

        # if exist_cap:
        #     lims = getzlims(
        #         pd.read_excel(fin, sheet_name=tech_type + "_exist_z", skiprows=0),
        #         techs,
        #         zones,
        #     )

        #     if lims.size != 0:
        #         param_outdd.append(
        #             getzlims(
        #                 pd.read_excel(
        #                     fin, sheet_name=tech_type + "_exist_z", skiprows=0
        #                 ),
        #                 techs,
        #                 zones,
        #             )
        #         )

        if esys_cap:
            dout = params.loc[
                (~params["Esys_capacity~2050"].isnull()), "Esys_capacity~2050"
            ]

            if dout.empty:
                continue
            else:
                wrapdd(
                    data2dd(
                        dout.values,
                        [
                            params[
                                (params["Technology Name (highRES)"] != "pgen")
                                & (~params["Esys_capacity~2050"].isnull())
                            ]["Technology Name (highRES)"].values
                        ],
                        rounddp=2,
                    ),
                    tech_type + "_fx_natcap",
                    root / out / (esys + "_" + tech_type + "_fx_natcap.dd"),
                )

        for s in sets:
            if s == "g" or s == "s":
                vals = params["Technology Name (highRES)"]

                set_outdd.append(wrapdd(data2dd(vals, []), s, "set"))
                continue

                # outsets.append(wrapdd(data2dd(vals,[])))
                # set_outdd.append(wrapdd(data2dd(vals,[],s,"set"))
            elif s == "non_vre" or s == "vre" or s == "hydro_res":
                vals = params[params["set"] == s]["Technology Name (highRES)"]

                set_outdd.append(
                    wrapdd(data2dd(vals, []), s + "(" + sets[0] + ")", "set")
                )
                
                
            else:
                vals = params.loc[(~params[s].isnull()), "Technology Name (highRES)"]

                set_outdd.append(
                    wrapdd(
                        data2dd(vals, []),
                        tech_type + "_" + s + "(" + sets[0] + ")",
                        "set",
                    )
                )
                
        for s in data_out:
            for p in data_out[s]:
                if "varom" in p or "fuel" in p:
                    # rounddp=5
                    rounddp = 8
                elif "emis" in p:
                    # rounddp=3
                    rounddp = 8
                else:
                    # rounddp=2
                    rounddp = 8

                out_par = tech_type + "_" + str.replace(p, " ", "")

                if "capex" in p or "fuelcost" in p:
                    out_par = re.sub(r"\d+", "", out_par)

                vals = params.loc[~(params[p].isnull()), p].values
                t_out = params.loc[
                    ~(params[p].isnull()), "Technology Name (highRES)"
                ].values

                param_outdd.append(
                    wrapdd(
                        data2dd(vals, [t_out], rounddp=rounddp), out_par, "parameter"
                    )
                )

        param_outdd = np.concatenate(param_outdd, axis=0)
        set_outdd = np.concatenate(set_outdd, axis=0)

        pad = np.repeat(np.array(""), set_outdd.shape[0]).reshape(set_outdd.shape[0], 1)

        outdd = np.concatenate((np.hstack((set_outdd, pad)), param_outdd), axis=0)

        np.savetxt(
            root / out / (run + "_" + tech_type + ".dd"), outdd, delimiter=" ", fmt="%s"
        )

        # param_outdd=np.concatenate(param_outdd,axis=0)
        # set_outdd=np.concatenate(set_outdd,axis=0)

        # np.savetxt(
        #     root/out/(run+"_"+tech_type+"_parameters.dd"), param_outdd,
        #     delimiter=" ",
        #     fmt="%s")
        # np.savetxt(
        #     root/out/(run+"_"+tech_type+"_sets.dd"),
        #     set_outdd, delimiter=" ",
        #     fmt="%s")


def euro_demand2dd(
    europedemandcsvlocation,
    aggregated_regions,
    dpath,
    opath,
    dstart,
    dstop,
    scen_db,
    esys_scen,
    yr,
):
    d = pd.read_csv(europedemandcsvlocation)

    try:
        d["datetime"] = pd.to_datetime(d["datetime"], format="%d/%m/%Y %H:%M")
    except ValueError:
        d["datetime"] = pd.to_datetime(d["datetime"], format="%Y-%m-%d %H:%M:%S")

    d = d.set_index("datetime")

    d = d.loc[:, d.columns.isin(aggregated_regions)]

    d = d[dstart:dstop]

    # TODO a better warning than below

    if d.shape[1] != len(aggregated_regions):
        print("Countries missing...")

    d[(d == 0)] = np.nan

    d.interpolate(limit=2, inplace=True)

    if pd.isnull(d).any().any():
        for c in d.columns[pd.isnull(d).any(axis=0)]:
            print(d.loc[pd.isnull(d[c]), c])

        print("Zero demands for: ", d.columns[pd.isnull(d).any(axis=0)])

    t = np.arange(d.shape[0])
    z = d.columns.values

    wrapdd(
        data2dd(d.values.T, [z, t], all_combin=True),
        "demand",
        "parameter",
        outfile=opath / (esys_scen + "_demand_" + str(yr) + ".dd"),
    )
