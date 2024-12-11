from allib import general as gl
from allib import hindtoolcalc as hc_calc
from allib import hindtoolplot as hc_plt
from allib import latex as ltx

# %% Startup Block
import os
import argparse
import inspect
import datetime
import shutil
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

script_name = os.path.basename(__file__)

parser = argparse.ArgumentParser(description=script_name)
parser.add_argument('-i', metavar='path_in', required=False, type=str, default='Input.txt',
                    help='the filepath to the input file, if empty "Input.txt"')
parser.add_argument('-o', metavar='path_out', required=False, type=str,
                    help='the filepath to the output dir, if empty, taken from Input file')

args = parser.parse_args()

path_in = args.i

filename = inspect.getframeinfo(inspect.currentframe()).filename
path_main = os.path.dirname(os.path.abspath(filename))

print(f"reading Inputfile ({path_in})...")

INPUT = gl.read_input_txt(path_in)
timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

if args.o is None:
    if INPUT['General']['dir_name'] is None:
        path_out = os.path.abspath(INPUT['General']['path_out']) + '\\HindCast_' + timestamp + '\\'
    else:
        path_out = os.path.abspath(INPUT['General']['path_out']) + '\\' + INPUT['General']['dir_name'] + '\\'

else:
    path_out = os.path.abspath(args.o) + '/'

print(f"Path_out = {path_out}")

if not os.path.exists(path_out):
    os.makedirs(path_out)

shutil.copy(path_in, path_out + 'Input.txt')

# %% program specific
# Database and output configurations
db_path = INPUT["General"]["db_path"]
table_name = 'Hind_combined'

columnnames = list(INPUT["ColumNames"].values())
# %% calculate

if INPUT['General']['calculate']:

    print("calculate...")
    beta_range = INPUT["Parameters"]["beta"]
    # alpha_range = INPUT["Parameters"]["alpha"]

    Data_Out = {}
    # for alpha in alpha_range:
    for beta in beta_range:

        print(f'   beta = {beta}')
        Calc = hc_calc.Calculation()
        df = Calc.initilize_from_db(db_path, table_name, columnnames)

        Calc.add_filter(mode='nans')
        Calc.apply_filters()

        T_p = df[INPUT["ColumNames"]["T_p"]]
        H_s = df[INPUT["ColumNames"]["H_s"]]
        dir_wave = df[INPUT["ColumNames"]["dir_T_mean"]]
        dir_wind = df[INPUT["ColumNames"]["dir_v_m"]]
        v_m = df[INPUT["ColumNames"]["v_m"]]

        indizes_swell, indizes_wind, rating = gl.separate_wind_swell(
            T_p, v_m, dir_wave, dir_wind, INPUT["Parameters"]["d"], None, 1, beta
        )

        result = {}

        result['indizes_swell'] = indizes_swell
        result['indizes_wind'] = indizes_wind
        result['rating'] = rating
        result['parameters'] = {'beta': beta}

        Seg = hc_calc.Segment(0, angles=None,
                              result=result,
                              angle_name=None,
                              colnames={'Hs': H_s.name, 'Tp': T_p.name, 'vm': v_m.name},
                              indizes=list(df.index))

        Calc.result = [Seg]

        Data_Out[f"beta={beta}"] = Calc

# %% write to DB
if INPUT["General"]["write_to_db"]:
    print("write to database...")

    for Calc_name, Calc in Data_Out.items():
        Seg = Calc.result[0]

        Tiles = []

        df = Calc.load_from_db(colnames_ini=True)
        T_p = df[INPUT["ColumNames"]["T_p"]]
        H_s = df[INPUT["ColumNames"]["H_s"]]
        v_m = df[INPUT["ColumNames"]["v_m"]]
        titels = Calc.create_segment_title()

        T_p_wind = T_p.loc[Seg.result["indizes_wind"]]
        T_p_wind.name = f"T_p_wind, WaveAge, beta={Seg.result['parameters']['beta']}"

        T_p_swell = T_p.loc[Seg.result["indizes_swell"]]
        T_p_swell.name = f"T_p_swell, WaveAge, beta={Seg.result['parameters']['beta']}"

        H_s_wind = H_s.loc[Seg.result["indizes_wind"]]
        H_s_wind.name = f"H_s_wind, WaveAge, beta={Seg.result['parameters']['beta']}"

        H_s_swell = H_s.loc[Seg.result["indizes_swell"]]
        H_s_swell.name = f"H_s_swell, WaveAge, beta={Seg.result['parameters']['beta']}"

        df = pd.concat([T_p_wind, T_p_swell, H_s_wind, H_s_swell], axis=1, join='outer')

        gl.add_dataframe_to_db(db_path, "Hind_combined", df)
        gl.add_dataframe_to_db(db_path, "Hind_raw_WindSwell_Seperation", df)

    META = gl.export_df_from_sql(db_path, "Hind_MetaData")
    META.loc["WindSwell_Seperation", "Source"] = "JBO intern preprocessing"
    META.loc["WindSwell_Seperation", "Date created"] = str(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    META.loc["WindSwell_Seperation", "Water Depth"] = INPUT["Parameters"]["d"]
    META.loc["WindSwell_Seperation", "Start Date"] = str(df.index[0])
    META.loc["WindSwell_Seperation", "End Date"] = str(df.index[-1])
    META.loc["WindSwell_Seperation", "Time Step"] = (df.index[1] - df.index[0]).total_seconds()
    META.loc["WindSwell_Seperation", "Number of samples"] = len(df)

    gl.add_dataframe_to_db(db_path, "Hind_MetaData", META)

# %% Plot
if INPUT["General"]["plot"]:
    print("plot...")
    figsize_fullpage = [size * 0.39370079 for size in INPUT["General"].get("writing_box", {})]
    figsize_fullpage_caption = [figsize_fullpage[0], figsize_fullpage[1] * 0.9]

    figsize_halfpage = [figsize_fullpage[0], figsize_fullpage[1] / 2.5]
    figsize_thirdpage = [figsize_fullpage[0], figsize_fullpage[1] / 3]
    figsize_twothirdpage = [figsize_fullpage[0], figsize_fullpage[1] / 1.5]
    figsize_halfpage_halfpage = [figsize_fullpage[0] / 2, figsize_fullpage[1] / 2.5]

    for Calc_name, Calc in Data_Out.items():
        print(f"    {Calc_name}")
        Seg = Calc.result[0]

        Tiles = []

        df = Calc.load_from_db(colnames_ini=True)
        T_p = df[INPUT["ColumNames"]["T_p"]]
        H_s = df[INPUT["ColumNames"]["H_s"]]
        v_m = df[INPUT["ColumNames"]["v_m"]]

        titels = Calc.create_segment_title()

        T_p_wind = T_p.loc[Seg.result["indizes_wind"]].values
        T_p_swell = T_p.loc[Seg.result["indizes_swell"]].values

        v_m_wind = v_m.loc[Seg.result["indizes_wind"]].values
        v_m_swell = v_m.loc[Seg.result["indizes_swell"]].values

        H_s_wind = H_s.loc[Seg.result["indizes_wind"]].values
        H_s_swell = H_s.loc[Seg.result["indizes_swell"]].values

        rating = Seg.result["rating"]

        # hstp wind
        tile_hstp_wind = hc_plt.Tile(num=1,
                                     x_label=INPUT["Aliase"]["H_s_wind"],
                                     y_label=INPUT["Aliase"]["T_p_wind"],
                                     title=f"$T_p(H_s)$ Wind sea, {round(len(H_s_wind) / len(H_s) * 100, 1)} \\%")

        scatter = hc_plt.Scatter(x=H_s_wind,
                                 y=T_p_wind,
                                 size=5,
                                 cmap='cool',
                                 cmap_norm='sqrt')

        tile_hstp_wind.add_scatter(scatter)
        Tiles.append(tile_hstp_wind)

        # vmHs wind
        tile_vmhs_wind = hc_plt.Tile(num=1,
                                     x_label=INPUT["Aliase"]["v_m"],
                                     y_label=INPUT["Aliase"]["H_s_wind"],
                                     title=f"$H_s(v_m)$ Wind sea, {round(len(H_s_wind) / len(H_s) * 100, 1)} \\%")

        scatter = hc_plt.Scatter(x=v_m_wind,
                                 y=H_s_wind,
                                 size=5,
                                 cmap='cool',
                                 cmap_norm='sqrt')

        tile_vmhs_wind.add_scatter(scatter)
        Tiles.append(tile_vmhs_wind)

        # hstp swell
        tile_hstp_swell = hc_plt.Tile(num=1,
                                      x_label=INPUT["Aliase"]["H_s_swell"],
                                     y_label=INPUT["Aliase"]["T_p_swell"],
                                     title=f"$T_p(H_s)$ Swell sea, {round(len(H_s_swell) / len(H_s) * 100, 1)} \\%")

        scatter = hc_plt.Scatter(x=H_s_swell,
                                 y=T_p_swell,
                                 size=5,
                                 cmap='cool',
                                 cmap_norm='sqrt')

        tile_hstp_swell.add_scatter(scatter)
        Tiles.append(tile_hstp_swell)

        # vmhs swell
        tile_vmhs_swell = hc_plt.Tile(num=1,
                                      x_label=INPUT["Aliase"]["v_m"],
                                     y_label=INPUT["Aliase"]["H_s_wind"],
                                     title=f"$H_s(v_m)$ Swell sea, {round(len(H_s_swell) / len(H_s) * 100, 1)} \\%")

        scatter = hc_plt.Scatter(x=v_m_swell,
                                 y=H_s_swell,
                                 size=5,
                                 cmap='cool',
                                 cmap_norm='sqrt')

        tile_vmhs_swell.add_scatter(scatter)
        Tiles.append(tile_vmhs_swell)

        FIG = hc_plt.plot_tiled(Tiles, global_max=[None, None], grid=[2,2], global_min=[0, 0], figsize=figsize_twothirdpage, use_pgf=False)
        gl.save_figs_as_png(FIG, path_out + Calc_name + '_wind_swell', dpi=300)

        Tiles = []
        # hstp total
        tile_hstp_total = hc_plt.Tile(num=1,
                                      x_label=INPUT["Aliase"]["H_s"],
                                      y_label=INPUT["Aliase"]["T_p"],
                                      title=f"$T_p(H_s)$ Total sea")

        scatter = hc_plt.Scatter(x=H_s.values,
                                 y=T_p.values,
                                 size=5,
                                 cmap='cool',
                                 cmap_norm='sqrt')

        tile_hstp_total.add_scatter(scatter)
        Tiles.append(tile_hstp_total)

        # vmhs total
        tile_vmhs_total = hc_plt.Tile(num=1,
                                      x_label=INPUT["Aliase"]["v_m"],
                                      y_label=INPUT["Aliase"]["H_s"],
                                      title=f"$H_s(v_m)$ Total sea")

        scatter = hc_plt.Scatter(x=v_m.values,
                                 y=H_s.values,
                                 size=5,
                                 cmap='cool',
                                 cmap_norm='sqrt')

        tile_vmhs_total.add_scatter(scatter)
        Tiles.append(tile_vmhs_total)

        # hstp rating
        tile_hstp_rating = hc_plt.Tile(num=1,
                                       x_label=INPUT["Aliase"]["H_s"],
                                      y_label=INPUT["Aliase"]["T_p"],
                                      title=f"$T_p(H_s)$ Total sea, wave age criterion $\\beta$",
                                       scatter_min=0)

        scatter = hc_plt.Scatter(x=H_s.values,
                                 y=T_p.values,
                                 size=5,
                                 cmap_mode='manual',
                                 c=rating,
                                 cmap='cool',
                                 alpha=0.5,
                                 cbar=True,
                                 cbar_label="$\\beta = \\frac{U_{10}}{c} \cos \left(\\theta-\\theta_w\\right)$")

        cbar_extraticks = [(Seg.result["parameters"]["beta"], 'beta')]
        tile_hstp_rating.add_scatter(scatter)
        Tiles.append(tile_hstp_rating)

        # vmhs rating
        tile_vmhs_rating = hc_plt.Tile(num=1,
                                       x_label=INPUT["Aliase"]["v_m"],
                                      y_label=INPUT["Aliase"]["H_s"],
                                      title=f"$H_s(v_m)$ Total sea, wave age criterion $\\beta$",
                                       scatter_min=0)

        scatter = hc_plt.Scatter(x=v_m.values,
                                 y=H_s.values,
                                 size=5,
                                 cmap_mode='manual',
                                 c=rating,
                                 cmap='cool',
                                 alpha=0.5,
                                 cbar=True,
                                 cbar_label="$\\beta=\\frac{U_{10}}{c} \cos \left(\\theta-\\theta_w\\right)$")

        tile_vmhs_rating.add_scatter(scatter)
        Tiles.append(tile_vmhs_rating)

        FIG = hc_plt.plot_tiled(Tiles, global_max=[None, None], global_min=[0, 0], grid=[2,2], scatter_max=None, scatter_min=None, figsize=figsize_twothirdpage, use_pgf=False)
        gl.save_figs_as_png(FIG, path_out + Calc_name + '_rating', dpi=300)




# %% Report
if INPUT["General"]["Report"]:

    # Create Report Output
    path_report = path_out + "report"
    try:
        # Create the new folder
        os.makedirs(path_report, exist_ok=True)  # exist_ok=True prevents an error if the folder already exists
        print(f"Folder created successfully at: {path_report}")
    except Exception as e:
        print(f"An error occurred: {e}")

    TEMPLATES = ltx.load_templates(os.path.abspath(".\\latex_templates"))

    # load figures in Dataframe (from output dir or optional dir)
    if INPUT["General"]["fig_path"] is not None:
        path_figs = os.path.abspath(INPUT["General"]["fig_path"])
    else:
        path_figs = os.path.abspath(path_out)

    FIGURES = ltx.load_figures(path_figs)

    TEX = {}
    # Data Basis
    chapter = 'WindSwellSeperation'

    TEX[chapter] = TEMPLATES[chapter]

    TEX[chapter] = ltx.insertLatexVars(TEX[chapter], {'beta': INPUT["Parameters"]["beta_report"],
                                                      'd': INPUT["Parameters"]["d"]})

    Fig_rating = FIGURES.loc[f"beta={INPUT['Parameters']['beta_report']}_rating_page_1"]
    Fig_rating["caption"] = "total sea density and wave age $\\beta$"

    Fig_wind_swell = FIGURES.loc[f"beta={INPUT['Parameters']['beta_report']}_wind_swell_page_1"]
    Fig_wind_swell["caption"] = "wind - swell seperation"

    TEX[chapter] = ltx.include_Fig(TEX[chapter], Fig_rating)
    TEX[chapter] = ltx.include_Fig(TEX[chapter], Fig_wind_swell)

    for name, tex in TEX.items():
        with open(path_report + r'\\' + name + '.tex', 'w', encoding='utf-8') as file:
            file.write(tex)