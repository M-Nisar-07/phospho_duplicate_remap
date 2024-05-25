import pandas as pd
df = pd.read_csv("duplicate_sites.csv")

print(df.shape)

df = df.groupby(["pmid","new_uniprot","mapped_genesymbol","exp_condition"]).agg(pd.Series.tolist)

print(df.shape)

result_dict = df.set_index('mapped_phosphosite')['expression'].to_dict()
print(result_dict)
df.to_excel("grouped.xlsx")