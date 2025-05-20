import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("mapping_coverage_HP126.txt", sep="\t", header=None, names=["contig", "pos", "depth"])
plt.plot(df["pos"], df["depth"])
plt.xlabel("Genomic Position")
plt.ylabel("Read Depth")
plt.title("Coverage Profile")
plt.show()
