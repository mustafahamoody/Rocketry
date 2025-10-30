import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("LSF5_data.csv")

# print(df.head())
plt.figure(1)
plt.plot(df["time"], df["[OPT-102] Ox Tank"])
plt.plot(df["time"], df["[CCPT-101] Combustion Chamber PT1 (psi)"])
plt.xlabel("Time (ms)")
plt.ylabel("Ox Tank Pressure (psi)")
plt.figure(2)
plt.plot(df["[CCPT-101] Combustion Chamber PT1 (psi)"]/df["[OPT-102] Ox Tank"])
plt.xlabel("Time (ms)")
plt.ylabel("CC %")
plt.show()