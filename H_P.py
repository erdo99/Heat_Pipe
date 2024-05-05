import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import pandas as pd
import numpy as np
import math
import csv
from tkinter import StringVar

import matplotlib.pyplot as plt
from scipy.integrate import quad

# Tkinter ile arayüzü burada oluşturuluyor
root = tk.Tk()
root.configure(bg="lightyellow")
p = root.title("Heat Pipe Calculator")
root.geometry("300x800")
pi = math.pi


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
class HP:

  def __init__(self):
    self.Material = {
        'WATER': {
            'Sigma_hat': 1.0,
            'M': 0.018,
            'k': 0.6,
            'Psat': {
                'T_data': [
                    273.1500, 275.1500, 277.1500, 283.1500, 287.1500, 291.1500,
                    293.1500, 298.1500, 303.1500, 307.1500, 313.1500, 317.1500,
                    323.1500, 327.1500, 333.1500, 343.1500, 353.1500, 363.1500,
                    369.1500, 373.1500, 383.1500, 393.1500, 403.1500, 413.1500,
                    423.1500
                ],
                'Pro_data': [
                    0.0061, 0.0071, 0.0081, 0.0123, 0.0160, 0.0206, 0.0234,
                    0.0317, 0.0425, 0.0533, 0.0738, 0.0911, 0.1235, 0.1502,
                    0.1995, 0.3120, 0.4741, 0.7018, 0.8777, 1.0142, 1.4338,
                    1.9867, 2.7028, 3.6154, 4.7616
                ]
            },
            'Sigma': {
                'Pro_data': [
                    0.0756,
                    0.0749,
                    0.0742,
                    0.0728,
                    0.0712,
                    0.0696,
                    0.0679,
                    0.0662,
                    0.0644,
                    0.0626,
                    0.0608,
                    0.0589,
                    0.0482,
                ],
                'T_data': [
                    273.1500,
                    278.1500,
                    283.1500,
                    293.1500,
                    303.1500,
                    313.1500,
                    323.1500,
                    333.1500,
                    343.1500,
                    353.1500,
                    363.1500,
                    373.1500,
                    423.1500,
                ]
            },
            'Vis': {
                'T_data': [
                    273.1500,
                    283.1500,
                    293.1500,
                    298.1500,
                    303.1500,
                    313.1500,
                    323.1500,
                    333.1500,
                    343.1500,
                    353.1500,
                    363.1500,
                    373.1500,
                    383.1500,
                    393.1500,
                    413.1500,
                    433.1500,
                ],
                'Pro_data': [
                    1.0e-05 * 0.1792,
                    1.0e-05 * 0.1306,
                    1.0e-05 * 0.1004,
                    1.0e-05 * 0.0893,
                    1.0e-05 * 0.0801,
                    1.0e-05 * 0.0658,
                    1.0e-05 * 0.0553,
                    1.0e-05 * 0.0474,
                    1.0e-05 * 0.0413,
                    1.0e-05 * 0.0364,
                    1.0e-05 * 0.0326,
                    1.0e-05 * 0.0294,
                    1.0e-05 * 0.0268,
                    1.0e-05 * 0.0246,
                    1.0e-05 * 0.0212,
                    1.0e-05 * 0.0188,
                ]
            },
            'hfg': {
                'T_data': [
                    273.1500,
                    275.1500,
                    277.1500,
                    283.1500,
                    287.1500,
                    291.1500,
                    293.1500,
                    298.1500,
                    303.1500,
                    307.1500,
                    313.1500,
                    317.1500,
                    323.1500,
                    327.1500,
                    333.1500,
                    343.1500,
                    353.1500,
                    363.1500,
                    369.1500,
                    373.1500,
                    383.1500,
                    393.1500,
                    413.1500,
                    433.1500,
                ],
                'Pro_data': [
                    2500900,
                    2496200,
                    2491400,
                    2477200,
                    2467700,
                    2458300,
                    2453500,
                    2441700,
                    2429800,
                    2420300,
                    2406000,
                    2396400,
                    2381900,
                    2372300,
                    2357700,
                    2333000,
                    2308000,
                    2282500,
                    2266900,
                    2256400,
                    2229600,
                    2202100,
                    2144300,
                    2082000,
                ]
            }
        },
        'IPA': {
            'Sigma_hat': 1.0,
            'M': 0.0601,
            'k': 0.14,
            'Psat': {
                'T_data': [
                    280,
                    282,
                    284,
                    286,
                    288,
                    290,
                    292,
                    294,
                    296,
                    298,
                    300,
                    302,
                    304,
                    306,
                    308,
                    310,
                    312,
                    314,
                    316,
                    318,
                    320,
                    322,
                    324,
                    326,
                    328,
                    330,
                    332,
                    334,
                    336,
                    338,
                    340,
                    342,
                    344,
                    346,
                    348,
                    350,
                    352,
                    354,
                    356,
                    358,
                    360,
                    362,
                    364,
                    366,
                    368,
                    370,
                    372,
                    374,
                    376,
                    378,
                    380,
                ],
                'Pro_data': [
                    1700, 2000, 2300, 2600, 3000, 3500, 3900, 4500, 5100, 5800,
                    6500, 7400, 8300, 9300, 10400, 11700, 13100, 14600, 16300,
                    18100, 20100, 22300, 24700, 27300, 30100, 33200, 36500,
                    40100, 44000, 48200, 52800, 57700, 63000, 68600, 74700,
                    81200, 88200, 95700, 103700, 112200, 121300, 131000,
                    141300, 152200, 163800, 176100, 189100, 202900, 217500,
                    232900, 249200
                ]
            },
            'Sigma': {
                'Pro_data': [
                    0.0229,
                    0.0223,
                    0.0217,
                    0.0212,
                    0.0207,
                    0.0202,
                    0.0197,
                    0.0192,
                    0.0187,
                    0.0182,
                    0.0174,
                    0.0172,
                    0.0169,
                    0.0165,
                    0.0162,
                    0.0159,
                ],
                'T_data': [
                    273.1500,
                    287.1500,
                    293.1500,
                    298.1500,
                    303.1500,
                    308.1500,
                    313.1500,
                    318.1500,
                    323.1500,
                    328.5500,
                    336.7500,
                    340.1500,
                    344.1500,
                    347.9500,
                    350.6500,
                    353.4500,
                ]
            },
            'Vis': {
                'T_data': [
                    283.1500, 293.1500, 298.1500, 303.1500, 313.1500, 323.1500,
                    335.1500, 344.3500, 354.6500
                ],
                'Pro_data': [
                    1.0e-05 * 0.4190, 1.0e-05 * 0.3080, 1.0e-05 * 0.2650,
                    1.0e-05 * 0.2310, 1.0e-05 * 0.1760, 1.0e-05 * 0.1380,
                    1.0e-05 * 0.1030, 1.0e-05 * 0.0875, 1.0e-05 * 0.0695
                ]
            },
            'hfg': {
                'T_data': [
                    280,
                    282,
                    284,
                    286,
                    288,
                    290,
                    292,
                    294,
                    296,
                    298,
                    300,
                    302,
                    304,
                    306,
                    308,
                    310,
                    312,
                    314,
                    316,
                    318,
                    320,
                    322,
                    324,
                    326,
                    328,
                    330,
                    332,
                    334,
                    336,
                    338,
                    340,
                    342,
                    344,
                    346,
                    348,
                    350,
                    352,
                    354,
                    356,
                    358,
                    360,
                    362,
                    364,
                    366,
                    368,
                    370,
                    372,
                    374,
                    376,
                    378,
                    380,
                ],
                'Pro_data': [
                    775990, 774020, 772010, 770000, 767870, 765780, 763540,
                    761340, 759100, 756790, 754440, 752020, 749550, 747130,
                    744550, 741910, 739320, 736570, 733760, 731000, 728090,
                    725120, 722190, 719110, 715870, 712580, 709340, 706040,
                    702690, 699290, 695740, 692130, 688580, 684880, 681120,
                    677220, 673380, 669480, 665440, 661350, 657220, 653050,
                    648830, 644570, 640160, 635820, 631330, 626810, 622250,
                    617640, 613000
                ]
            }
        },
        'METHANOL': {
            'Sigma_hat': 1.0,
            'M': 0.032,
            'k': 0.2010,
            'Psat': {
                'T_data': [
                    263.1500, 268.1500, 273.1500, 278.1500, 283.1500, 288.1500,
                    293.1500, 298.1500, 303.1500, 308.1500, 313.1500, 318.1500,
                    323.1500, 328.1500, 333.1500, 338.1500, 343.1500, 348.1500,
                    353.1500, 358.1500, 363.1500, 368.1500, 373.1500, 378.1500,
                    383.1500, 388.1500, 393.1500
                ],
                'Pro_data': [
                    2019, 2844, 3947, 5404, 7304, 9755, 12881, 16826, 21757,
                    27864, 35362, 44493, 55527, 68763, 84531, 103194, 125146,
                    150816, 180667, 215199, 254947, 300483, 352417, 411397,
                    478109, 553279, 637674
                ]
            },
            'Sigma': {
                'Pro_data': [
                    0.0263, 0.0236, 0.0226, 0.0225, 0.0218, 0.0211, 0.0205,
                    0.0201, 0.0185, 0.0166, 0.0146, 0.0125, 0.0104
                ],
                'T_data': [
                    263.1500, 283.1500, 293.1500, 298.1500, 303.1500, 313.1500,
                    318.1500, 323.1500, 343.1500, 363.1500, 383.1500, 403.1500,
                    423.1500
                ]
            },
            'Vis': {
                'T_data': [
                    248.2000, 273.2000, 283.1500, 293.1500, 298.1500, 303.1500,
                    308.1500, 313.1500, 318.1500, 323.1500, 348.2000, 373.2000
                ],
                'Pro_data': [
                    1.0e-05 * 0.1508, 1.0e-05 * 0.0982, 1.0e-05 * 0.0839,
                    1.0e-05 * 0.0750, 1.0e-05 * 0.0699, 1.0e-05 * 0.0650,
                    1.0e-05 * 0.0620, 1.0e-05 * 0.0575, 1.0e-05 * 0.0554,
                    1.0e-05 * 0.0516, 1.0e-05 * 0.0400, 1.0e-05 * 0.0321
                ]
            },
            'hfg': {
                'T_data': [
                    223.1500, 243.1500, 263.1500, 283.1500, 289.1500, 303.1500,
                    323.1500, 343.1500, 363.1500, 383.1500, 403.1500, 423.1500
                ],
                'Pro_data': [
                    1194000, 1187000, 1182000, 1175000, 1167000, 1155000,
                    1125000, 1085000, 1035000, 980000, 920000, 850000
                ]
            }
        }
    }
    self.Geometry = None
    self.Solver = {
        'Tinf': None,
        'hconv': None,
        'alpha': None,
        'Tv': None,
        'Qpc': 0.2
    }
    self.Results = {
        'z': [],
        'R': [],
        'Theta': [],
        'R_final': [],
        'q_axial': [],
        'T_final': [],
        'm_dot_final': [],
        'Theta_final': [],
        'Tevap': None,
        'Tcon': None,
        'Qc': None,
        'Reff': None,
        'keff': None,
        'Tv': None,
        'q_pc': 0.2
    }


Empty_HP = HP()

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#labela = tk.Label(root, text="Heat input : ")
#labela.pack(pady=20, padx=20)
frame = tk.Frame(root, bg="lightyellow")
#my_entry1 = tk.Entry(frame, width=10)
#my_entry1.grid(row=0, column=1, sticky=tk.W + tk.E)

#labelb = tk.Label(root, text="Wick type :", font=("Arial", 10))
#labelb.pack(padx=2, pady=2)

#my_entry2 = tk.Entry(frame, width=20)
#my_entry2.grid(row=1, column=1, sticky=tk.W + tk.E)

frame.columnconfigure(0, weight=1)
frame.columnconfigure(1, weight=1)
frame.columnconfigure(2, weight=2)
frame.columnconfigure(3, weight=2)
frame.columnconfigure(4, weight=2)
frame.columnconfigure(5, weight=2)
frame.columnconfigure(6, weight=2)
frame.columnconfigure(7, weight=2)
frame.columnconfigure(8, weight=2)
frame.columnconfigure(9, weight=2)
frame.columnconfigure(10, weight=2)
label1 = tk.Label(frame, text="Heat input : ", anchor=tk.W)
label1.grid(row=0, column=0, sticky=tk.W + tk.E)
label2 = tk.Label(frame, text="Coolant Temperature : ", anchor=tk.W)
label2.grid(row=1, column=0, sticky=tk.W + tk.E)
label3 = tk.Label(frame, text="Groove height : ", anchor=tk.W)
label3.grid(row=2, column=0, sticky=tk.W + tk.E)

label4 = tk.Label(frame, text="Groove width : ", anchor=tk.W)
label4.grid(row=3, column=0, sticky=tk.W + tk.E)

label5 = tk.Label(frame, text="Fin width : ", anchor=tk.W)
label5.grid(row=4, column=0, sticky=tk.W + tk.E)

label6 = tk.Label(frame, text="Number of grooves : ", anchor=tk.W)
label6.grid(row=5, column=0, sticky=tk.W + tk.E)

label7 = tk.Label(frame, text="Heat sink length : ", anchor=tk.W)
label7.grid(row=6, column=0, sticky=tk.W + tk.E)

label8 = tk.Label(frame, text="Heat Source length : ", anchor=tk.W)
label8.grid(row=7, column=0, sticky=tk.W + tk.E)

#label9 = tk.Label(frame, text="Heat Source length : ")
#label9.grid(row=8, column=0, sticky=tk.W + tk.E)

label10 = tk.Label(frame, text="Total length : ", anchor=tk.W)
label10.grid(row=8, column=0, sticky=tk.W + tk.E)

label11 = tk.Label(frame, text="Base Thickness : ", anchor=tk.W)
label11.grid(row=9, column=0, sticky=tk.W + tk.E)

label12 = tk.Label(frame, text="Heat Sink : ", anchor=tk.W)
label12.grid(row=10, column=0, sticky=tk.W + tk.E)

label13 = tk.Label(frame, text="Contact Angle : ", anchor=tk.W)
label13.grid(row=11, column=0, sticky=tk.W + tk.E)

label14 = tk.Label(frame, text="Thermal Conductivity : ", anchor=tk.W)
label14.grid(row=12, column=0, sticky=tk.W + tk.E)

label15 = tk.Label(frame, text="Working Fluid : ", anchor=tk.W)
label15.grid(row=13, column=0, sticky=tk.W + tk.E)

label16 = tk.Label(frame, text="Wick Type : ", anchor=tk.W)
label16.grid(row=14, column=0, sticky=tk.W + tk.E)
entries = []


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
def Show_Value():
  global T_sat_new
  T_sat_new = 0.2
  global q_pc_new
  q_pc_new = 0.2
  Geometry = []
  for entry in entries:
    if entry.get() == "":
      messagebox.showerror("Error", "Please enter all the values")
      return None
  # groove width along heat pipe [m]
  global W_g
  W_g = float(entries[3].get()) * 1e-3
  # groove/fin height [m]
  global H_g
  H_g = float(entries[2].get()) * 1e-3
  # fin width along heat pipe [m]
  global W_f
  W_f = float(entries[4].get()) * 1e-3
  # Number of grooves in HP
  global N_g
  N_g = float(entries[5].get())
  #print(str(N_g))
  # condenser length[m]
  global L_c
  L_c = float(entries[6].get()) * 1e-3

  # evaporator length[m]
  global L_e
  L_e = float(entries[7].get()) * 1e-3

  # Total length of HP [m]
  global L_t
  L_t = float(entries[8].get()) * 1e-3
  #print(str(L_t))
  # base height [m]
  global H_b
  H_b = float(entries[9].get()) * 1e-3
  # width of the side area after grooves [m]
  global W_s
  W_s = 0.5 * 1e-3

  # adiabatic length[m]
  global L_a
  L_a = L_t - L_c - L_e
  # Coolant Temperature [K]
  global CoolTemp
  CoolTemp = float(entries[1].get()) + 273.15
  # heat sink convection coefficient [W/m^2K]
  global hconv
  hconv = float(entries[10].get())
  # heat input  [W]
  global Qtotal
  Qtotal = float(entries[0].get())
  # Thermal Conductivity [W/mK]
  global ThermalCond
  ThermalCond = float(entries[12].get())
  # Contact Angle [deg]
  global Alpha
  Alpha = float(entries[11].get()) * pi / 180
  #print(str(Alpha))
  Empty_HP.Solver['Qpc'] = float(0.8 * Qtotal)
  #print(str(Empty_HP.Solver['Qpc']))

  selected_fluid = selected_fluid_var.get()
  if selected_fluid == "WATER":
    HP_Material = Empty_HP.Material["WATER"]
  elif selected_fluid == "METHANOL":
    HP_Material = Empty_HP.Material["METHANOL"]
  elif selected_fluid == "IPA":
    HP_Material = Empty_HP.Material["IPA"]
  else:
    messagebox.showerror("Error", "Seçilen Sıvı Tanımlı Değil!")
  Empty_HP.Solver['Tv'] = float(CoolTemp + 5)  #****!!!!!!!!!!
  #print(str(Empty_HP.Solver['Tv']))
  HPT()  #HESAPLAMA FONKSİYONU BURADA ÇAĞRILIYOR

  color1 = [0, 0.4470, 0.7410]
  color2 = [0.8500, 0.3250, 0.0980]

  # İlk veri seti
  plt.plot(Empty_HP.Results['z'] * 1e3,
           Empty_HP.Results['R'] * 1e3,
           color=color1,
           linewidth=2)
  plt.xlabel('z [mm]', fontname='Times New Roman', fontsize=20)
  plt.ylabel(r'$\mathcal{R}\;\rm{[mm]}$',
             fontname='Times New Roman',
             fontsize=20)
  plt.gca().tick_params(axis='y')
  plt.gca().yaxis.set_major_formatter('{:.2f}'.format)

  # İkinci veri seti
  plt.twinx()
  plt.plot(Empty_HP.Results['z'][1 + 10::10] * 1e3,
           Empty_HP.Results['T_final'] - 273.15,
           color=color2,
           linewidth=2)
  plt.ylabel('Wall temperature [Celcius]',
             fontname='Times New Roman',
             fontsize=20)

  # Diğer GUI bileşenlerine değer atama
  Tmax = max(Empty_HP.Results['T_final']) - 273.15
  Tmin = min(Empty_HP.Results['T_final']) - 273.15
  qpc_value = Empty_HP.Results['Qpc']
  qc_value = Empty_HP.Results['Qc']
  Reff_value = Empty_HP.Results['Reff']
  keff_value = Empty_HP.Results['keff'] / 1e3

  plt.show()


#?????????????????????????????????????????????????????????????
button = tk.Button(root, text="Show", command=Show_Value, bg="light blue")

button.pack(padx=10, pady=20, side=tk.BOTTOM)

for i in range(13):
  my_entry = tk.Entry(frame, width=20)

  my_entry.grid(row=i, column=1, sticky=tk.W + tk.E)

  entries.append(my_entry)
selected_fluid_var = StringVar()
comboBox2 = ttk.Combobox(frame, values=["Grooved", "Sintered"])
comboBox2.grid(row=14, column=1, sticky=tk.W + tk.E)
comboBox2.current(0)
comboBox1 = ttk.Combobox(frame,
                         textvariable=selected_fluid_var,
                         values=["IPA", "METHANOL", "WATER"])
comboBox1.grid(row=13, column=1, sticky=tk.W + tk.E)
comboBox1.current(0)


# Heat Pipe ın ana hesaplama fonksiyonu
def HPT():
  R_u = 8.314  #Universal gas constant [J/mol.K]
  selected_fluid = selected_fluid_var.get()
  HP_Material = Empty_HP.Material[selected_fluid]

  M = HP_Material["M"]
  #print(str(M))
  K = HP_Material["k"]
  Ks = ThermalCond
  #print(str(Ks))
  Sigma_hat = HP_Material["Sigma_hat"]
  #print(str(Sigma_hat))
  Corr_tot = 1.056 - 1.056 * ((pi - 2 * Alpha) / 2 / math.cos(Alpha)**2 - math.tan(Alpha)) / \
  (((pi - 2 * Alpha) / math.cos(Alpha)**2 - 2 * math.tan(Alpha)) + 8 * H_g / W_g - (pi - 2 * Alpha) / math.cos(Alpha)**2 + 2 * math.tan(Alpha))
  #print(str(Corr_tot))
  #*********************************************************************
  # Sliced Properties
  dl = 0.001
  #print(str(dl))
  n = math.floor(L_e / dl)  #number of slice in evaporator
  #print(str(n) + " bu n dir.")
  #print(str(L_e / dl))
  m = math.floor(L_a / dl)
  #print(str(m) + " Bu m dir.")  # number of slice in adibatic
  p = math.floor(L_c / dl)
  Nslice = n + m + p
  res = 10  # resolution in every slice
  Ltotal = L_e + L_a + L_c
  # **********************************************
  T_sat = Empty_HP.Solver['Tv']
  Q_pc = Empty_HP.Solver['Qpc']
  Q = Qtotal  # Total Heat input [W]
  Q_in = np.ones(n) * Q / n  # Heat input per slice [W]
  Q_pc_e = np.ones(
      n
  ) * Q_pc / n  # Estimated phase change heat transfer per slice in evaporator [W]
  Q_pc_a = np.zeros(
      m)  # Estimated phase change heat transfer per slice in adibatic [W]
  Q_pc_c = -np.ones(
      p
  ) * Q_pc / p  # Estimated phase change heat transfer per slice in condenser [W]
  q_pc = np.concatenate((Q_pc_e, Q_pc_a, Q_pc_c))
  # print(str(q_pc))
  T_inf = CoolTemp  # Ambientte Temperature [K]
  h_conv = hconv  # Heat transfer coefficent at heat sink [W/m^2K]
  tol = 1e-6
  tol_CO = 1e-10
  error_q_pc = 1.0
  error_T_sat = 1.0
  iter = 1
  epsf = 0.3  # should be 0 < epsf < 1
  epsg = 0.9  # should be 0 < epsg < 1
  N_int = 500  # number of division for integration at thin film
  N_int_g = 500  # number of division for integration at groove film
  N_int_ft = 500  # number of division for integration at fin top
  z = np.arange(0, Ltotal + Ltotal / (res * Nslice), Ltotal / (res * Nslice))

  # print(str(z))
  z_e = z[0:(res * n + 1)]
  z_a = z[(res * n + 1):(res * (n + m) + 1)]
  z_c = z[(res * (n + m) + 1):(res * (n + m + p) + 1)]
  while error_T_sat > tol or error_q_pc > tol:
    h_fg = np.interp(T_sat, HP_Material['hfg']['T_data'],
                     HP_Material['hfg']['Pro_data'])
    P_sat = np.interp(T_sat, HP_Material['Psat']['T_data'],
                      HP_Material['Psat']['Pro_data'])
    Sigma = np.interp(T_sat, HP_Material['Sigma']['T_data'],
                      HP_Material["Sigma"]['Pro_data'])
    vis = np.interp(T_sat, HP_Material['Vis']['T_data'],
                    HP_Material['Vis']['Pro_data'])

    # Calculation of Coeffiecient a from Kibetic Theory
    a = 2 * Sigma_hat / (2 - Sigma_hat) * np.sqrt(M / (2* np.pi * R_u * T_sat)) * \
    (M * P_sat * h_fg) / (R_u * T_sat * (T_sat +1))
    Kappa = np.pi**6 * vis * (W_g**2 + 4 * H_g**2) / (256 * W_g**3 * H_g**3 *
                                                      Corr_tot)
    #***************************************************************************
    #Total m_dot for each slice
    m_dot = np.zeros(n + m + p)
    for i in range(1, n + m + p):
      m_dot[i - 1] = np.sum(q_pc[:i] / (N_g * h_fg))

#***********************************************
# Evaporator region
    E = np.zeros(3 * n - 1)
    #print(str(len(E)))
    E[0] = n * m_dot[0] / L_e
    E[1] = P_sat - Sigma / (W_g / (2 * np.cos(Alpha)))
    E[2] = n * (m_dot[1] - m_dot[0]) / L_e
    E[3] = m_dot[0]
    E[4] = Kappa * E[0] * L_e**2 / (n**2 * 2) + E[1]
    E[5::3] = n * (m_dot[2:n] - m_dot[1:n - 1]) / L_e
    E[6::3] = m_dot[1:n - 1]
    for i in range(3, 3 * n - 4, 3):
      E[i + 4] = Kappa * E[i - 1] * L_e**2 / (
          n**2 * 2) + Kappa * E[i] * L_e / n + E[i + 1]
    P_e = np.zeros(res * n + 1)
    m_dot_e = np.zeros(res * n + 1)
    P_e[0:res + 1] = Kappa * E[0] * z_e[0:res + 1] ** 2 / 2 + E[1]
    m_dot_e[0:res + 1] = E[0] * z_e[0:res + 1]

    for j in range(1, n):
        tmp = res * (j - 1) + 1
        P_e[tmp:(tmp + res)] = (Kappa * E[3 * j - 3] * (z_e[tmp:(tmp + res)] - (j - 1) * L_e / n) ** 2) / 2 + \
                               Kappa * E[3 * j - 2] * (z_e[tmp:(tmp + res)] - (j - 1) * L_e / n) + E[3 * j - 1]
        m_dot_e[tmp:(tmp + res)] = E[3 * j - 3] * (z[tmp:(tmp + res)] - (j - 1) * L_e / n) + E[3 * j - 2]

    Theta_e = np.arccos((P_sat - P_e) * W_g / (2 * Sigma))

    print("Üçüncü P_e bu: ")
    print(P_e)
    #print(Theta_e)
    print("Kappa : ")
    print(Kappa)
    print(" z_e : ")
    print(z_e)
    print("Sigma : ")
    print(Sigma)

    Radiee_e = Sigma / (P_sat - P_e)
    Theta_e_av = np.zeros(n)
    for i in range(1, n + 1):
      start_index = res * (i - 1)
      end_index = res * i
      Theta_e_av[i - 1] = np.sum(Theta_e[start_index:end_index]) / res

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Adiabatic Region
    D = np.zeros(3 * m)
    print(str(len(D)))
    D[0] = m * (m_dot[n] - m_dot[n - 1]) / L_a
    D[1] = m_dot[n]
    D[2] = Kappa * E[3 * n - 4] * L_e**2 / (n**2 * 2) + Kappa * E[
        3 * n - 3] * L_e / n + E[3 * n - 2]  # değiştirildi
    D[3::3] = m * (m_dot[n + 1:n + m] - m_dot[n:n + m - 1]) / L_a
    P_a = np.zeros(res * m)
    m_dot_a = np.zeros(res * m)

    for j in range(1, m + 1):
      tmp = res * (j - 1)
      P_a[tmp:tmp + res] = (Kappa * D[3 * j - 3] * (z_a[tmp:tmp + res] - L_e - \
                                                    (j - 1) * L_a / m)**2) / 2 + \
      Kappa * D[3 * j -2] * (z_a[tmp: tmp + res] - L_e - (j-1) * L_a / m) + D[3*j-1]
      m_dot_a[tmp: tmp + res] = D[3 * j -2] * (z_a[tmp:tmp + res] - L_e - (j-1) *
      L_a/m) + \
      D[3*j -1]     # değiştirildi index 144 is out of bounds for axis 0 with size
    Theta_a = np.arccos((P_sat - P_a) * W_g / (2 * Sigma))
    Radii_a = Sigma / (P_sat - P_a)
    Theta_a_av = np.zeros(m)
    for i in range(1, m + 1):
      start_index = res * (i - 1)
      end_index = res * i
      Theta_a_av[i - 1] = np.sum(Theta_a[start_index:end_index])
#???????????????????????????????????????????????????
# Condenser Region
    C = np.zeros(3 * p)
    C[0] = p * (m_dot[n + m] - m_dot[n + m - 1]) / L_c
    C[1] = m_dot[n + m]
    C[2] = Kappa * D[3 * m - 2] * L_a**2 / (m**2 * 2) + Kappa * D[3 * m -
                                                                  1] * L_a / m
    C[3::3] = p * (m_dot[n + m + 1:n + m + p] -
                   m_dot[n + m:n + m + p - 1]) / L_c
    C[4::3] = m_dot[n + m:n + m + p - 1]
    for i in range(1, 3 * p - 3, 3):
      C[i + 4] = Kappa * C[i - 1] * L_c**2 / (
          p**2 * 2) + Kappa * C[i] * L_c / p + C[i + 1]
    P_c = np.zeros(res * p)
    m_dot_c = np.zeros(res * p)
    for j in range(1, p + 1):
      tmp = res * (j - 1)
      P_c[tmp:tmp + res] = (Kappa * C[3*j - 3] * (z_c[tmp:tmp + res] - L_e - L_a - (j-1) * L_c/p)**2) / 2 + \
       Kappa * C[3*j - 2] * (z_c[tmp:tmp + res] - L_e - L_a - (j-1) * L_c/p) + \
       C[3*j-1]
      m_dot_c[tmp:tmp + res] = C[3*j -2] *(z_c[tmp:tmp + res] - L_e - L_a - (j-1) *
      L_c/p) + \
      C[3 * j -1 ]
    Theta_c = np.arccos((P_sat - P_c) * W_g / (2 * Sigma))
    Radii_c = Sigma / (P_sat - P_c)
    Theta_c_av = np.zeros(p)
    for i in range(1, p + 1):
      start_index = res * (i - 1)
      end_index = res * i
      Theta_c_av[i - 1] = np.sum(Theta_c[start_index:end_index]) / res
      #**********************************************************************
    P = np.concatenate((P_e, P_a, P_c))
    Theta = np.concatenate((Theta_e, Theta_a, Theta_c))
    Radii = np.concatenate((Radiee_e, Radii_a, Radii_c))

    A_side = 2 * W_s * (H_g + H_b)
    A_base = N_g * (W_f * H_g + (W_f + W_g) * H_b) + A_side
    R_a_ax = np.ones(m + 1)
    R_e = np.zeros(n)

    # axial base resistance [K/Wm]
    R_e_ax = np.zeros(n)
    R_e_ax[:n - 1] = L_e / (n * Ks * A_base)
    R_a_ax[0] = (L_e / (2 * n) + L_a / (2 * m)) / (Ks * A_base)
    R_a_ax[1:m] = L_a / (m * Ks * A_base)
    R_a_ax[-1] = (L_c / (2 * p) + L_a / (2 * m)) / (Ks * A_base)
    R_c_ax = np.zeros(p - 1)
    R_c_ax[:] = L_c / (p * Ks * A_base)

    R_ax = np.concatenate((R_e_ax, R_a_ax, R_c_ax))
    #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    #EVAPORATOR RESISTANCE CALCULATION
    R = 0.5 * W_g / np.cos(Theta_e_av)
    # R_tf^e (fin tarafı sıvı ve faz  değişim direncinin hesaplanması)
    x0 = R * np.sin(Theta_e_av)
    xend = epsf * (R - x0) + x0
    CC0_evaporator_store = np.zeros(n)
    for i in range(n):
      x = np.linspace(x0[i], xend[i], N_int)
      gamma = np.arcsin(x / R[i]) - Theta_e_av[i]
      S = R[i] * gamma
      I = W_g / 2 * (x - x0[i] ) -0.5 * R[i] **2 * (np.arcsin(x / R[i]) - \
                                                  np.arcsin(x0[i] / R[i] )) \
                    -0.5 * ( x * np.sqrt(R[i] **2 - x**2) - \
                            x0[i] * np.sqrt(R[i] **2 - \
                                           x0[i] **2))
      delta_e = (I[1:] - I[:-1]) / (x[1:] - x[:-1])
      R_tf_l_e = delta_e / (K * (x[1:] - x[:-1]) * L_e / n)
      R_tf_pc_e = 1 / (a * h_fg * (S[:-1]) * L_e / n)
      R_tf_e = np.sum(1 / (R_tf_l_e + R_tf_pc_e))**(-1)
      #Calculation of R_g(Groove side liquid and phase change resistance)
      y0 = 0
      yend = R[i] * np.sin(np.arccos(xend[i] / R[i]))
      y = np.linspace(y0, yend, N_int_g)
      beta = np.arcsin(y / R[i])
      S = R[i] * beta
      h1 = epsg * (H_g - R[i] + x0[i])
      I = H_g * (y - y0) + x0[i] * (y - y0) -h1 * (y - y0)  \
          -0.5 * R[i] **2 * (np.arcsin(y / R[i]) - np.arcsin(y0 / R[i])) - \
                           0.5 * (y * np.sqrt(R[i] **2 -y**2) \
                                  -y0 * np.sqrt(R[i] **2 - y0 **2))

      delta_e_g = (I[1:] - I[:-1]) / (y[1:] - y[:-1])
      R_g_l_e = delta_e_g / (K * (y[1:] - y[:-1]) * L_e / n)
      R_g_pc_e = 1 / (a * h_fg * (S[1:] - S[:-1]) * L_e / n)
      R_g_e_top = np.sum(1 / (R_g_l_e + R_g_pc_e))**(-1)
      #*************************************************************************
      R_fb_e = (H_b + H_g) / (Ks * L_e / n * 0.5 * W_f
                              )  # Base resistance fin side [K/W]
      R_gb_e = H_b / (Ks * L_e / n * 0.5 * W_g
                      )  # Base resistance fin side [K/W]
      R_g_l_e_base = h1 / (K * L_e / n * 0.5 * W_g
                           )  # Liquid resistance groove side [K/W]

      R_g_e = R_gb_e + R_g_l_e_base + R_g_e_top  # Groove resistance [K/W]
      if Q_pc_e[i] >= 0:
        R_f_e = R_fb_e + R_tf_e  # Fin resistance [K/W]
        R_e[i] = (1 / R_f_e + 1 / R_g_e)**(-1) / (2 * N_g)
      else:
        # *************************************************
        # Calculation of R_ft^c
        # *************************************************
        # 4th-order polynomial for the film profile on the top of the fin
        CC0 = 0.5e-6
        CC1 = -np.tan(0.5 * np.pi - Theta_e_av[i])
        CC2 = 0
        CC3 = -2 * CC1 / W_f**2
        CC4 = -4 * CC1 / W_f**3
        error_C0 = 1

        def error_mdotC(CC0):

          def int_fun(s):
            return (a * (-0.01)) / (1 + (a *
                                         (CC0 + CC1 * (s - 0.5 * W_f) + CC2 *
                                          (s - 0.5 * W_f)**2 + CC3 *
                                          (s - 0.5 * W_f)**3 + CC4 *
                                          (s - 0.5 * W_f)**4) * h_fg) / K)

          integral_result = quad(int_fun, 0, 0.5 * W_f)[0]
          return Sigma / (3 * vis) * 6 * (CC0**3) * CC3 + integral_result

        while abs(error_C0) > tol_CO:
          error_C0 = error_mdotC(CC0)
          CC0 = CC0 - 0.01 * CC0 * error_mdotC(CC0) / (
              error_mdotC(CC0 + 0.01 * CC0) - error_mdotC(CC0))

          pass
        CC0_evaporator_store[i] = CC0

        w_f0 = 0
        w_fend = 0.5 * W_f
        w = np.linspace(w_f0, w_fend, N_int_ft)
        I_ft = CC0 * (w - w_f0) + CC1 / 2 * (w - 0.5 * W_f) ** 2 - CC1 / 2 * (w_f0 - 0.5 * W_f) ** 2 + \
        CC2 / 3 * (w - 0.5 * W_f) ** 3 - CC2 / 3 * (w_f0 - 0.5 * W_f) ** 3 + \
        CC3 / 4 * (w - 0.5 * W_f) ** 4 - CC3 / 4 * (w_f0 - 0.5 * W_f) ** 4 + \
        CC4 / 5 * (w - 0.5 * W_f) ** 5 - CC4 / 5 * (w_f0 - 0.5 * W_f) ** 5

        delta_ft = np.diff(I_ft) / np.diff(w)

        R_ft_l_e = delta_ft / (K * np.diff(w) * L_e / n)

        def fun(e):
          return np.sqrt(1 + (CC1 + 2 * CC2 * (e - 0.5 * W_f) + 3 * CC3 *
                              (e - 0.5 * W_f)**2 + 4 * CC4 *
                              (e - 0.5 * W_f)**3)**2)

        sfintop, _ = quad(fun, 0, 0.5 * W_f)
        R_ft_pc_e = 1 / (a * h_fg * (sfintop / N_int_ft) * L_e / n)

        # R_ft_e hesaplama
        R_ft_e = np.sum(1 / (R_ft_l_e + R_ft_pc_e))**(-1)

        # R_f_e hesaplama
        R_f_e = R_fb_e + (1 / R_tf_e + 1 / R_ft_e)**(-1)

        # R_e hesaplama
        R_e[i] = (1 / R_f_e + 1 / R_g_e)**(-1) / (2 * N_g)

      pass
    # ADIABATIC REGION RESISTANCE CALCULATION
    R_a = np.zeros(m)
    R = 0.5 * W_g / np.cos(Theta_a_av)

    # Calculation of R_tf^e (Fin side liquid and phase change resistance)
    x0 = R * np.sin(Theta_a_av)
    xend = epsf * (R - x0) + x0
    CC0_adiabatic_store = np.zeros(m)
    for i in range(m):
      x = np.linspace(x0[i], xend[i], N_int)
      gamma = np.arcsin(x / R[i]) - Theta_a_av[i]
      S = R[i] * gamma

      I = W_g / 2 * (x - x0[i]) - 0.5 * R[i]**2 * (np.arcsin(x / R[i]) - np.arcsin(x0[i] / R[i])) - \
          0.5 * (x * np.sqrt(R[i]**2 - x**2) - x0[i] * np.sqrt(R[i]**2 - x0[i]**2))

      delta_a = (I[1:] - I[:-1]) / (x[1:] - x[:-1])

      R_tf_l_a = delta_a / (K * (x[1:] - x[:-1]) * L_a / m)
      R_tf_pc_a = 1 / (a * h_fg * (S[1:] - S[:-1]) * L_a / m)

      R_tf_a = np.sum(1 / (R_tf_l_a + R_tf_pc_a))**(-1)
      # *************************************************
      # Calculation of R_g (Groove side liquid and phase change resistance)
      y0 = 0
      yend = R[i] * np.sin(np.arccos(xend[i] / R[i]))
      y = np.linspace(y0, yend, N_int_g)
      beta = np.arcsin(y / R[i])
      S = R[i] * beta
      h1 = epsg * (H_g - R[i] + x0[i])
      I = H_g * (y - y0) + x0[i] * (y - y0) - h1 * (y - y0) - 0.5 * R[i]**2 * (np.arcsin(y / R[i]) - np.arcsin(y0 / R[i])) - \
          0.5 * (y * np.sqrt(R[i]**2 - y**2) - y0 * np.sqrt(R[i]**2 - y0**2))
      delta_a_g = np.diff(I) / np.diff(y)
      R_g_l_a = delta_a_g / (K * np.diff(y) * L_a / m)
      R_g_pc_a = 1 / (a * h_fg * np.diff(S) * L_a / m)
      R_g_a_top = np.sum(1 / (R_g_l_a + R_g_pc_a))**(-1)
      # Baz dirençlerini hesaplama
      R_fb_a = (H_b + H_g) / (Ks * (L_a / m) *
                              (0.5 * W_f))  # Taban direnci fin tarafı [K/W]
      R_gb_a = H_b / (Ks * (L_a / m) *
                      (0.5 * W_g))  # Taban direnci yiv tarafı [K/W]
      R_g_l_a_base = h1 / (K * (L_a / m) * 0.5 * W_g
                           )  # Sıvı direnci yiv tarafı [K/W]

      R_g_a = R_gb_a + R_g_l_a_base + R_g_a_top  # Yiv direnci [K/W]
      if Q_pc_a[i] >= 0:
        # Fin direncini hesaplama
        R_f_a = R_fb_a + R_tf_a  # Fin direnci [K/W]

        R_a[i] = 1 / (1 / R_f_a + 1 / R_g_a) / (2 * N_g)

      else:
        # Fin üstü film profilinin 4. dereceden polinom ile hesaplanması
        CC0 = 5e-6
        CC1 = -np.tan(0.5 * np.pi - Theta_a_av[i])
        CC2 = 0
        CC3 = -2 * CC1 / W_f**2
        CC4 = -4 * CC1 / W_f**3
        error_C0 = 1

        def int_fun(s):
          return (a * (-0.01)) / (1 + (a * (CC0 + CC1 * (s - 0.5 * W_f) + CC2 *
                                            (s - 0.5 * W_f)**2 + CC3 *
                                            (s - 0.5 * W_f)**3 + CC4 *
                                            (s - 0.5 * W_f)**4) * h_fg) / K)

        # Hata fonksiyonu
        def error_mdotC(CC0):
          integral_result, _ = quad(int_fun, 0, 0.5 * W_f)
          return (Sigma / (3 * vis)) * 6 * (CC0**3) * CC3 + integral_result

        while abs(error_C0) > tol_CO:
          current_error = error_mdotC(CC0)
          derivative_error = (error_mdotC(CC0 + 0.01 * CC0) -
                              current_error) / (0.01 * CC0)
          error_C0 = current_error
          CC0 -= 0.01 * CC0 * error_C0 / derivative_error
        CC0_adiabatic_store[i, 0] = CC0

        w_f0 = 0
        w_fend = 0.5 * W_f
        w = np.linspace(w_f0, w_fend, N_int_ft)

        I_ft = CC0 * (w - w_f0) + \
               CC1 / 2 * (w - 0.5 * W_f) ** 2 - CC1 / 2 * (w_f0 - 0.5 * W_f) ** 2 + \
               CC2 / 3 * (w - 0.5 * W_f) ** 3 - CC2 / 3 * (w_f0 - 0.5 * W_f) ** 3 + \
               CC3 / 4 * (w - 0.5 * W_f) ** 4 - CC3 / 4 * (w_f0 - 0.5 * W_f) ** 4 + \
               CC4 / 5 * (w - 0.5 * W_f) ** 5 - CC4 / 5 * (w_f0 - 0.5 * W_f) ** 5
        delta_ft = np.diff(I_ft) / np.diff(w)

        # R_ft_l_a hesaplama
        R_ft_l_a = delta_ft / (K * np.diff(w) * L_a / m)

        # fun ve sfintop hesaplama
        def fun(e):
          return np.sqrt(1 + (CC1 + 2 * CC2 * (e - 0.5 * W_f) + 3 * CC3 *
                              (e - 0.5 * W_f)**2 + 4 * CC4 *
                              (e - 0.5 * W_f)**3)**2)

        sfintop, _ = quad(fun, 0, 0.5 * W_f)

        # R_ft_pc_a hesaplama
        R_ft_pc_a = 1 / (a * h_fg * (sfintop / N_int_ft) * L_a / m)

        # R_ft_a hesaplama
        R_ft_a = sum(1 / (R_ft_l_a + R_ft_pc_a))**(-1)
        R_f_a = R_fb_a + (1 / R_tf_a + 1 / R_ft_a)**(-1
                                                     )  # Fin resistance [K/W]

        R_a[i] = (1 / R_f_a + 1 / R_g_a)**(-1) / (2 * N_g)  # Geri dön düzelt
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #CONDENSER RESISTANCE CALCULATION
    R_c = np.zeros(p)
    R = 0.5 * W_g / np.cos(Theta_c_av)
    #Calculation of R_tf^e (Fin side liquid and phase change resistance)
    x0 = R * np.sin(Theta_c_av)
    xend = epsf * (R - x0) + x0
    CC0_condenser_store = np.zeros(p)
    for i in range(p):
      x = np.linspace(x0[i], xend[i], N_int)
      gamma = np.arcsin(x / R[i]) - Theta_c_av[i]
      S = R[i] * gamma

      I = W_g / 2 * (x - x0[i]) - 0.5 * R[i]**2 * (np.arcsin(x / R[i]) - np.arcsin(x0[i] / R[i])) - \
          0.5 * (x * np.sqrt(R[i]**2 - x**2) - x0[i] * np.sqrt(R[i]**2 - x0[i]**2))

      delta_c = (I[1:] - I[:-1]) / (x[1:] - x[:-1])

      R_tf_l_c = delta_c / (K * (x[1:] - x[:-1]) * L_c / p)
      R_tf_pc_c = 1 / (a * h_fg * (S[1:] - S[:-1]) * L_c / p)

      R_tf_c = np.sum(1 / (R_tf_l_c + R_tf_pc_c))**(-1)
      #Calculation of R_g (Groove side liquid and phase change resistance)
      y0 = 0
      yend = R[i] * np.sin(np.arccos(xend[i] / R[i]))
      y = np.linspace(y0, yend, N_int_g)
      beta = np.arcsin(y / R[i])
      S = R[i] * beta
      h1 = epsg * (H_g - R[i] + x0[i])

      I = (H_g) * (y - y0) + (x0[i]) * (y - y0) - h1 * (y - y0) - 0.5 * R[i]**2 * (np.arcsin(y / R[i]) - np.arcsin(y0 / R[i])) - \
          0.5 * (y * np.sqrt(R[i]**2 - y**2) - y0 * np.sqrt(R[i]**2 - y0**2))

      delta_c_g = (I[1:] - I[:-1]) / (y[1:] - y[:-1])

      R_g_l_c = delta_c_g / (K * (y[1:] - y[:-1]) * (L_c / p))
      R_g_pc_c = 1 / (a * h_fg * (S[1:] - S[:-1]) * (L_c / p))

      R_g_c_top = np.sum(1 / (R_g_l_c + R_g_pc_c))**(-1)
      #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      R_fb_c = (H_b + H_g) / (Ks * (L_c / p) *
                              (0.5 * W_f))  # Base resistance fin side [K/W]
      R_gb_c = H_b / (Ks * (L_c / p) *
                      (0.5 * W_g))  # Base resistance fin side [K/W]
      R_g_l_c_base = h1 / (K * (L_c / p) * 0.5 * W_g
                           )  # Liquid resistance groove side [K/W]

      R_g_c = R_gb_c + R_g_l_c_base + R_g_c_top  # Groove resistance [K/W]
      if Q_pc_c[i] >= 0:
        R_f_c = R_fb_c + R_tf_c  # Fin resistance [K/W]
        R_c[i] = (1 / R_f_c + 1 / R_g_c)**(-1) / (2 * N_g)
      else:
        # *************************************************
        # Calculation of R_ft^c
        # *************************************************
        # 4th order polynomial for the film profile on the top of the fin
        CC0 = 5e-6
        CC1 = -math.tan(0.5 * math.pi - Theta_c_av[i])
        CC2 = 0
        CC3 = -2 * CC1 / W_f**2
        CC4 = -4 * CC1 / W_f**3
        error_C0 = 1

        def int_fun(s):
          return (a * (-0.01)) / (1 + (a * (CC0 + CC1 * (s - 0.5 * W_f) + CC2 *
                                            (s - 0.5 * W_f)**2 + CC3 *
                                            (s - 0.5 * W_f)**3 + CC4 *
                                            (s - 0.5 * W_f)**4) * h_fg) / K)

        def error_mdotC(CC0):
          integrand, _ = quad(int_fun, 0, 0.5 * W_f)
          return (Sigma /
                  (3 * vis)) * 6 * (CC0**3) * (-2 * CC1 / W_f**2) + integrand

        while abs(error_C0) > tol_CO:
          error_C0 = error_mdotC(CC0)
          CC0 = CC0 - 0.01 * CC0 * error_mdotC(CC0) / (
              error_mdotC(CC0 + 0.01 * CC0) - error_mdotC(CC0))
        CC0_condenser_store[i] = CC0  # dik hale gelmesi gerekebilir.
        w_f0 = 0
        w_fend = 0.5 * W_f

        # w vektörünü oluşturun
        w = np.linspace(w_f0, w_fend, N_int_ft)

        # I_ft değerlerini hesaplayın
        I_ft = (CC0 * (w - w_f0) + CC1 / 2 * (w - 0.5 * W_f)**2 - CC1 / 2 *
                (w_f0 - 0.5 * W_f)**2 + CC2 / 3 * (w - 0.5 * W_f)**3 -
                CC2 / 3 * (w_f0 - 0.5 * W_f)**3 + CC3 / 4 *
                (w - 0.5 * W_f)**4 - CC3 / 4 * (w_f0 - 0.5 * W_f)**4 +
                CC4 / 5 * (w - 0.5 * W_f)**5 - CC4 / 5 * (w_f0 - 0.5 * W_f)**5)
        delta_ft = np.diff(I_ft) / np.diff(w)

        # R_ft_l_c değerlerini hesaplayın
        R_ft_l_c = delta_ft / (K * np.diff(w) * (L_c / p))

        def fun(e):
          return np.sqrt(1 + (CC1 + 2 * CC2 * (e - 0.5 * W_f) + 3 * CC3 *
                              (e - 0.5 * W_f)**2 + 4 * CC4 *
                              (e - 0.5 * W_f)**3)**2)

        sfintop, _ = quad(fun, 0, 0.5 * W_f)
        R_ft_pc_c = 1 / (a * h_fg * (sfintop / N_int_ft) * (L_c / p))

        R_ft_c = sum(1 / (R_ft_l_c + R_ft_pc_c))**(-1)

        # Fin resistance [K/W]
        R_f_c = R_fb_c + (1 / R_tf_c + 1 / R_ft_c)**(-1)

        # Overall resistance [K/W]
        R_c[i] = (1 / R_f_c + 1 / R_g_c)**(-1) / (2 * N_g)
    # Convective Heat Transfer Resistance
    Rconv = np.ones(p)
    Rconv_value = 1 / (h_conv * (N_g * W_g +
                                 (N_g + 1) * W_f + 2 * W_s) * L_c / p)
    Rconv *= Rconv_value
    #!!!!!!!!!!!!!!!!!!!!!!!!!!

    Q_in_part = [-Q_in]

    # SOLUTION OF THE LINEAR SYSTEM
    # Generation of the Matrices
    Q_in_array = Q_in_part[0]  # Q_in_part içindeki diziye erişim

    # Şimdi concatenate işlemini yapın
    A = np.concatenate(
        (Q_in_array, np.zeros(2 * n + 3 * m + 3 * p - 1), [-T_inf] * p, [0]))
    print("A Matrisi:", A)
    print(A.shape)

    # B Matrisi
    B_size = 3 * n + 3 * m + 4 * p  # B matrisinin boyutu
    B = np.zeros((B_size, B_size))
    print("B Matris boyutu  : ", B_size)
    print("B Matrisi :  ", B)

    # B1 - Birim matris (eye) B matrisinin bir bölümüne yerleştirilecek.
    B1_size = n + m + p
    B1 = -np.eye(B1_size - 1)
    B[0:B1_size - 1, 0:B1_size -
      1] = B1  # burasına -1 sonradan eklendi hatalı mı denemek lazım
    print("B Matrisi :  ", B)
    # B2 matrisi (burada kullanımı gösterilmemiş, dolayısıyla sadece tanımlıyorum)
    B2_size = n + m + p - 1

    B2 = np.zeros((n + m + p, n + m + p - 1))
    for i in range(n + m + p - 1):
      B2[i, i] = -1
      B2[i + 1, i] = 1
    print("B2 :", B2)
    B[:n + m + p, n + m + p:2 * n + 2 * m + 2 * p - 1] = B2
    B3 = np.eye(p)  # Birim matrisi oluşturduk

    B[n + m:n + m + p, 2 * n + 2 * m + 2 * p:2 * n + 2 * m + 3 * p] = B3
    B4 = -np.diag(R_ax)  # B4 matrisini oluşturduk  !! değiştirildi

    B[n + m + p:2 * n + 2 * m + 2 * p, n + m + p:2 * n + 2 * m + 2 * p] = B4
    B5_size = n + m + p - 1
    B5 = np.zeros((B5_size, n + m + p))

    # Değer atama döngüsü
    i = 0
    for j in range(B5_size):
      B5[i, j] = 1
      B5[i, j + 1] = -1
      i += 1
    B[n + m + p:2 * n + 2 * m + 2 * p - 1,
      2 * n + 2 * m + 3 * p:3 * n + 3 * m + 4 * p] = B5
    # B6 matrisi
    R_combined = np.concatenate((R_e, R_a, R_c))  # değiştirildi

    # -np.eye ile çarpma
    B6 = -np.eye(n + m + p) * R_combined[:, np.newaxis]
    B[2 * n + 2 * m + 2 * p:3 * n + 3 * m + 3 * p, :n + m + p] = B6

    # B7 matrisi
    B7 = np.eye(n + m + p)
    B[2 * n + 2 * m + 2 * p:3 * n + 3 * m + 3 * p,
      2 * n + 2 * m + 3 * p:3 * n + 3 * m + 4 * p] = B7

    # B8 matrisi
    B8 = -np.eye(
        p) * Rconv[:, np.
                   newaxis]  # Rconv vektörü için de benzer bir işlem yapılır.
    B[3 * n + 3 * m + 3 * p:3 * n + 3 * m + 4 * p,
      2 * n + 2 * m + 2 * p:2 * n + 2 * m + 3 * p] = B8

    # B9 matrisi
    B9 = -np.eye(p)
    B[3 * n + 3 * m + 3 * p:3 * n + 3 * m + 4 * p,
      3 * n + 3 * m + 3 * p:3 * n + 3 * m + 4 * p] = B9

    # B10 vektörü
    B10 = -np.ones((n + m + p, 1))
    B[2 * n + 2 * m + 2 * p:3 * n + 3 * m + 3 * p,
      3 * n + 3 * m + 4 * p - 1] = B10.squeeze()

    # B11 vektörü
    B11 = np.ones((1, n + m + p))
    B[3 * n + 3 * m + 4 * p - 1, 0:n + m + p] = B11

    print("B Matrisi : ")
    B_df = pd.DataFrame(B)
    B_df.to_csv('B_verileri.csv', index=False, header=True)
    #Solution of unknown vector x
    X_solve = np.linalg.solve(B, A.T)
    X = np.linalg.inv(B).dot(A.T)

    T_sat_new = X[3 * n + 3 * m + 4 * p]

    # q_pc_new vektörünü hesaplayalım ve parçalara ayıralım
    q_pc_new = X[:n + m + p].T
    q_pc_e = q_pc_new[:n]
    q_pc_a = q_pc_new[n:n + m]
    q_pc_c = q_pc_new[n + m:n + m + p]

    iter = iter + 1  # iterasyon sayısını artır

    # error_T_sat hesaplayalım
    error_T_sat = abs((T_sat_new - T_sat) / T_sat_new)

    # error_q_pc hesaplayalım
    error_q_pc = np.max(np.abs((q_pc_new - q_pc) / q_pc_new))
    if error_T_sat > tol or error_q_pc > tol:
      T_sat = T_sat_new
      q_pc = q_pc_new
      pass
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Empty_HP.Results['Tv'] = T_sat_new
  Empty_HP.Results['q_pc'] = q_pc_new
  Empty_HP.Results['z'] = z
  Empty_HP.Results['R'] = Radii
  Empty_HP.Results['Theta'] = Theta

  Empty_HP.Results['R_final'] = np.array([R_e, R_a, R_c]).T
  Empty_HP.Results['q_axial'] = X[n + m + p:2 * n + 2 * m + 2 * p - 1]
  Empty_HP.Results['T_final'] = X[2 * n + 2 * m + 3 * p:3 * n + 3 * m + 4 * p -
                                  1]
  Empty_HP.Results['m_dot_final'] = np.array([m_dot_e, m_dot_a, m_dot_c])
  Empty_HP.Results['Theta_final'] = np.array(Theta * 180 / np.pi).T
  Empty_HP.Results['Tevap'] = Empty_HP.Results['T_final'][0]
  Empty_HP.Results['Tcon'] = Empty_HP.Results['T_final'][-1]
  Empty_HP.Results['q_pc'] = np.sum(
      Empty_HP.Results['q_pc'] *
      (Empty_HP.Results['q_pc'] > 0))  # kontrol et akşam
  Empty_HP.Results['Qc'] = Qtotal - sum(
      Empty_HP.Results['q_pc'])  # dikkat et !!!!!!!!!
  Empty_HP.Results['Reff'] = (max(Empty_HP.Results['T_final']) -
                              min(Empty_HP.Results['T_final'])) / Qtotal
  Aeff = (Empty_HP.Results['Geometry']['w_g'] + Empty_HP.Results['Geometry']['w_f']) * Empty_HP.Results['Geometry'] ['N_g'] * \
   (Empty_HP.Results['Geometry']['h_g'] + Empty_HP.Results['Geometry']['h_b'])
  Empty_HP.Results[
      'keff'] = Empty_HP.Solver.Qtotal * Empty_HP.Geometry.Ltotal / Aeff / (
          max(Empty_HP.Results['T_final']) - min(Empty_HP.Results['T_final'])
      )  # sonradan dönmen gerekebilir


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

entries[0].insert(0, 13.2)
entries[1].insert(0, 25)
entries[2].insert(0, 0.8)
entries[3].insert(0, 0.8)
entries[4].insert(0, 0.8)
entries[0].focus()

entries[5].insert(0, 15)
entries[6].insert(0, 26)
entries[7].insert(0, 26)
entries[8].insert(0, 100)
entries[9].insert(0, 1.2)
entries[10].insert(0, 4000)
entries[11].insert(0, 15)
entries[12].insert(0, 180)

frame.pack(fill="x")
frame.columnconfigure(0, minsize=140)
frame.columnconfigure(1, minsize=30)

tk.mainloop()