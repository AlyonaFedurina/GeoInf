import requests
import math as m
import numpy as np
import matplotlib.pyplot as plt
from pyorbital.orbital import Orbital
from datetime import datetime, timedelta

#Переводим в радианы
def deg_to_rad(value):
    return (value * np.pi) / 180

#Переводим в градусы
def rad_to_deg(value):
    return (value * 180) / np.pi

#Переводим из сферических в декартовы координаты
def sph_to_dec(r, theta, phi):
    return r * m.cos(theta) * m.cos(phi), r * m.cos(theta) * m.sin(phi), r * m.sin(theta)

#Функция, задающая вектор посредством его координат и длины
def set_v(x1, y1, z1, x2, y2, z2):
    v = [x2 - x1, y2 - y1, z2 - z1]
    v_length = m.sqrt(V[0] ** 2 + V[1] ** 2 + V[2] ** 2)
    return v, v_length

#Считаем угол между векторами
def angle_v1_v2(v1, v2, v1_length, v2_length):
    return m.acos((v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) / (v1_length * v2_length))

#Считаем скалярное произведение векторов
def sc_pr(v1, v2):
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]

 #Получаем долготу, широту, высоту над Землёй спутника
def get_stll_data(TLE_1, TLE_2, utc_time):
    orb = Orbital("N", line1=TLE_1, line2=TLE_2)
    lo, lt, he = orb.get_lonlatalt(utc_time)
    return lo, la, he

#Получаем TLE
def get_TLE(file_with_TLE, stll_name):
    record = requests.get(file_with_TLE, stream=True)
    open('TLE.txt', 'wb').write(record.text)
    file = open('TLE.txt', 'r')
    temporary = file.read().split("\n")[:-1]
    for i in range(len(temporary)):

        if temporary[i] == stll_name:
            return [s tll_name, temporary[i + 1], temporary[i + 2]]

#Посчитаем расстояние от точки вне сферы до касательной
def dist_to_plane(x1, y1, z1, x2, y2, z2):
    dist = (x1 * x2 + y1 * y2 + z1 * z2 - x1 ** 2 - y1 ** 2 - z1 ** 2 ) / ((x1 ** 2 + y1 ** 2 + z1 ** 2) ** 0.5)
    return dist



he_LK = 0.197 #высота ЛК над морем
R = 6378.1375 #радиус Земли
la_LK_R = deg_to_rad(55.930096) #широта ЛК
lo_LK_R = deg_to_rad(37.517872) #долгота ЛК
x_LK, y_LK, z_LK = sph_to_dec(R + he_LK, la_LK_R,
                              lo_LK_R) #координаты ЛК в декартовой системе

D = (x_LK ** 2 + y_LK ** 2 + z_LK ** 2) #свободный член в уравнении касательной плоскости
z_pl = D / z_LK
y_pl = D / y_LK
x_pl = D / x_LK
North_V, North_V_L = set_v(0, 0, z_pl, x_LK, y_LK, z_LK)

Normal_V, Normal_V_L = set_v(0, 0, 0, x_LK, y_LK, z_LK)
East_V, East_V_L = set_v(0, 0, 0, North_V[1] * Normal_V[2] - Normal_V[1] * North_V[2],
                         -(North_V[0] * Normal_V[2] - North_V[2] * Normal_V[0]),
                         North_V[0] * Normal_V[1] - Normal_V[0] * North_V[1])

#Списки
arr_x = []
arr_y = []
arr_z = []
arr_time = []
arr_elevation = []
arr_azimuth = []
arr_usual_elevation = []
arr_usual_azimuth = []
arr_usual_time = []
tmp_elevtion_data = []
tmp_azimuth_data = []
tmp_time_data = []


list_TLE =  get_TLE("https://celestrak.com/NORAD/elements/active.txt", "NOAA 19                 ")
print(list_TLE)

print("Enter start date: MM HH DD MM YYYY: ")
start_time = input().split()
start_time = [int(c) for c in start_time]
start_time = [int(c) for c in start_time]
start_time = datetime(start_time[4], start_time[3], start_time[2], start_time[1], start_time[0])


print("Enter final date: MM HH DD MM YYYY: ")
end_time = input().split()
end_time = [int(c) for c in end_time]
end_time = datetime(end_time[4], end_time[3], end_time[2], end_time[1], end_time[0])
minutes = int((end_time - start_time).total_seconds() / 60)

for i in range(minutes):
    longitude, latitude, height = get_stll_data(list_TLE[1], list_TLE[2], start_time - timedelta(hours=3)) #перевод времени в международное
    start_time = start_time + timedelta(minutes=1)
    latitude_radians = deg_to_rad(latitude)
    longitude_radians = deg_to_rad(longitude)
    x, y, z = sph_to_dec(height + R, latitude_radians, longitude_radians)
    array_x.append(x)
    array_y.append(y)
    array_z.append(z)
    Sat_V, Sat_V_L = set_v(x_LK, y_LK, z_LK, x, y, z)
    dist_to_P = dist_to_plane(x_LK, y_LK, z_LK, x, y, z)
    dist_to_Sat = ((x - x_LK) ** 2 + (y - y_LK) ** 2 + (z - z_LK) ** 2) ** 0.5
    elevation = m.asin(dist_to_P / dist_to_Sat)
    azimuth = 0
    if elevation >= 0:
        Sat_V_N = [(Normal_V[0] * sc_pr(Normal_V, Sat_V)) / (Normal_V_L ** 2),
                   (Normal_V[1] * sc_pr(Normal_V, Sat_V)) / (Normal_V_L ** 2),
                   (Normal_V[2] * sc_pr(Normal_V, Sat_V)) / (Normal_V_L ** 2)]
        Sat_V_Pr, Sat_V_Pr_L = set_v(0, 0, 0, -Sat_V[0] + Sat_V_N[0], -Sat_V[1] + Sat_V_N[1], -Sat_V[2] + Sat_V_N[2])
        azimuth = angle_v1_v2(Sat_V_Pr, North_V, North_V_L, Sat_V_Pr_L)
        if 3 * np.pi / 2 >= angle_v1_v2(East_V, Sat_V_Pr, East_V_L, Sat_V_Pr_L) > np.pi / 2:
            azimuth = 2 * np.pi - azimuth
    elevation = rad_to_deg(elevation)
    array_azimuth.append(azimuth)

    array_elevation.append(elevation)

for i in range(len(array_time)):
    if 180 >= array_elevation[i] >= 0:
        temporary_elevation_data.append(array_elevation[i])
        temporary_azimuth_data.append(array_azimuth[i])
        temporary_time_data.append(array_time[i])
    else:
        if len(temporary_elevation_data) != 0:
            array_usual_time.append(temporary_time_data)
            array_usual_azimuth.append(temporary_azimuth_data)
            array_usual_elevation.append(temporary_elevation_data)
        temporary_elevation_data = []
        temporary_azimuth_data = []
        temporary_time_data = []

for i in range(len(array_usual_time)):
    print("[ Time: ", array_usual_time[i][0], "] [ Azimuth: ", rad_to_deg(array_usual_azimuth[i][0]), "] [ Elevation: ",
          max(array_usual_elevation[i]), "]")

print(array_usual_elevation)
print(array_usual_azimuth)

fig = plt.figure()
ax = fig.add_subplot(111, projection='polar')
ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)
ax.set_rlim(bottom=90, top=0)
for phi, theta in zip(array_usual_azimuth, array_usual_elevation):
    ax.plot(phi, theta)
fig.set_size_inches(7, 7)
plt.show()

sf = plt.figure()
ax = sf.add_subplot(111, projection='3d')
ax.plot(array_x, array_y, array_z)
ax.scatter(x_LK, y_LK, z_LK, color='black')
sf.set_size_inches(7, 7)
plt.show()
