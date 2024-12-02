from libaries import general as gl
import numpy as np
import matplotlib.pyplot as plt

# Database and output configurations
db_path = "C:/Users/aaron.lange/Desktop/Projekte/HindcastArchive/4_54.96112_6.254773_NIRAS/54.96112_6.254773_NIRAS.db"
table_name = 'Hind_combined'
column_names = ["Wind speed 10m LAT [m/s]", "Wind direction  10m LAT [degN-CF]", "Hs [m]", "Tp [s]", "Mean Wave Direction [degN-CF]"]
path_out = "C:/Users/aaron.lange/Desktop/Projekte/Hindcast_Tool/WindSwell_Seperator/test_3/"

# Load data
df = gl.export_df_from_sql(db_path, table_name, column_names=column_names)
df = df.dropna(how='any')

# Function to separate wind swell
def separate_wind_swell(T_p, v_m, dir_wave, dir_wind, water_depth, h_vm, alpha, beta):

    omega = 2 * np.pi / T_p
    k = gl.k_aus_omega(omega, water_depth)
    c = omega / k
    indizes_swell = []
    indizes_wind = []

    v_m = h_vm * v_m

    for T_p_curr, v_m_curr, dir_wave_curr, dir_wind_curr, c_curr, index in zip(
        T_p.values, v_m.values, dir_wave.values, dir_wind.values, c.values, T_p.index
    ):
        dir_wave_curr = dir_wave_curr * 2 * np.pi / 360
        dir_wind_curr = dir_wind_curr * 2 * np.pi / 360
        beta_compare = v_m_curr / c_curr * (np.cos(dir_wave_curr - dir_wind_curr))**alpha

        if beta_compare < beta:
            indizes_swell.append(index)
        else:
            indizes_wind.append(index)

    return indizes_swell, indizes_wind

# Parameters for alpha and beta
#alpha_range = np.linspace(0.1, 2, 20)  # Adjust range and number of elements
beta_range = np.linspace(0.6, 1.3, 10)
alpha = 2
# Data columns
T_p = df["Tp [s]"]
H_s = df["Hs [m]"]
dir_wave = df["Mean Wave Direction [degN-CF]"]
dir_wind = df["Wind direction  10m LAT [degN-CF]"]
v_m = df["Wind speed 10m LAT [m/s]"]

# Loop over alpha and beta combinations
#for alpha in alpha_range:
for beta in beta_range:
    print(f"{alpha}, {beta}")
    indizes_swell, indizes_wind = separate_wind_swell(
        T_p, v_m, dir_wave, dir_wind, 50, 150, alpha, beta
    )

    T_p_wind = T_p.loc[indizes_wind]
    T_p_swell = T_p.loc[indizes_swell]

    v_m_wind = v_m.loc[indizes_wind]
    v_m_swell = v_m.loc[indizes_swell]

    H_s_wind = H_s.loc[indizes_wind]
    H_s_swell = H_s.loc[indizes_swell]

    fig = plt.figure(figsize=(10, 12))
    fig.suptitle("sepration of wind swell" + "\n" + f"with alpha {alpha} and beta {beta}", fontsize=16)

    ax1 = fig.add_subplot(3, 2, 1)
    c = gl.c_scatterplot(H_s_wind, T_p_wind)
    ax1.scatter(H_s_wind, T_p_wind, c=v_m_wind, s=2)
    ax1.set_xlabel(H_s_wind.name)
    ax1.set_ylabel(T_p_wind.name)
    ax1.set_title(f"wind sea portion of points: {round(len(indizes_wind) / len(T_p) * 100)}")
    ax1.set_xlim(min(H_s), max(H_s))
    ax1.set_ylim(min(T_p), max(T_p))

    ax2 = fig.add_subplot(3, 2, 3)
    c = gl.c_scatterplot(H_s_swell, T_p_swell)
    ax2.scatter(H_s_swell, T_p_swell, c=v_m_swell, s=2)
    ax2.set_xlabel(H_s_swell.name)
    ax2.set_ylabel(T_p_swell.name)
    ax2.set_title(f"swell sea portion of points: {round(len(indizes_swell) / len(T_p) * 100)}")
    ax2.set_xlim(min(H_s), max(H_s))
    ax2.set_ylim(min(T_p), max(T_p))

    ax3 = fig.add_subplot(3, 2, 5)
    c = gl.c_scatterplot(H_s, T_p)
    ax3.scatter(H_s, T_p, c=v_m, s=2)
    ax3.set_xlabel(H_s.name)
    ax3.set_ylabel(T_p.name)
    ax3.set_title(f"total")
    ax3.set_xlim(min(H_s), max(H_s))
    ax3.set_ylim(min(T_p), max(T_p))

    # VMHS
    ax4 = fig.add_subplot(3, 2, 2)
    c = gl.c_scatterplot(v_m_wind, H_s_wind)
    ax4.scatter(v_m_wind, H_s_wind, c=c, s=2)
    ax4.set_xlabel(v_m_wind.name)
    ax4.set_ylabel(H_s_wind.name)
    ax4.set_title(f"wind sea portion of points: {round(len(indizes_wind) / len(T_p) * 100)}")
    ax4.set_xlim(min(v_m), max(v_m))
    ax4.set_ylim(min(H_s), max(H_s))

    ax5 = fig.add_subplot(3, 2, 4)
    c = gl.c_scatterplot(v_m_swell, H_s_swell)
    ax5.scatter(v_m_swell, H_s_swell, c=c, s=2)
    ax5.set_xlabel(v_m_swell.name)
    ax5.set_ylabel(H_s_swell.name)
    ax5.set_title(f"swell sea portion of points: {round(len(indizes_swell) / len(T_p) * 100)}")
    ax5.set_xlim(min(v_m), max(v_m))
    ax5.set_ylim(min(H_s), max(H_s))

    ax6 = fig.add_subplot(3, 2, 6)
    c = gl.c_scatterplot(v_m, H_s)
    ax6.scatter(v_m, H_s, c=c, s=2)
    ax6.set_xlabel(v_m.name)
    ax6.set_ylabel(H_s.name)
    ax6.set_title(f"total sea")
    ax6.set_xlim(min(v_m), max(v_m))
    ax6.set_ylim(min(H_s), max(H_s))

    plt.tight_layout()

    gl.save_figs_as_png([fig], path_out + f'alpha={alpha:.2f}_beta={beta:.2f}', dpi=300)


