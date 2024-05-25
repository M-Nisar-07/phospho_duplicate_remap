import pandas as pd
import numpy as np
def do(row):
    mapped_phosphosites = row['mapped_phosphosite']
    expressions = row['expression']
    map_dict = {}
    for phosphosite, expression in zip(mapped_phosphosites, expressions):
        if phosphosite not in map_dict:
            map_dict[phosphosite] = []
        map_dict[phosphosite].append(expression)
    return map_dict

def doo(x):
	li = []
	for val in x.values():
		if len(set(val)) == 1:
			li.append('True')
		else:
			li.append('False')
	return li


df = pd.read_csv("duplicate_sites.csv")


final_cols = list(df.columns)
FINAL_dF = pd.DataFrame(columns=final_cols)

grouped_df = df.groupby(["pmid", "new_uniprot", "mapped_genesymbol", "exp_condition"]).agg(
    lambda x: list(x)
).reset_index()

grouped_df['mappedval'] = grouped_df.apply(do, axis=1)
grouped_df['Boolean'] = grouped_df['mappedval'].apply(doo)

grouped_df["filter"] = grouped_df["Boolean"].apply(lambda x:True if "False" in x else False)


# ===========================Nisar code=============================================================



df = grouped_df.loc[grouped_df["filter"] == True ]



def getindex(x):
    value_indices = {}
    
    # Group the indices of identical values
    for index, item in enumerate(x):
        if item not in value_indices:
            value_indices[item] = []
        value_indices[item].append(index)
    
    # Convert the dictionary values to a list of lists
    separated_indices = list(value_indices.values())
    
    return separated_indices


df["indexes"] = df["mapped_phosphosite"].apply(getindex)

df = df.explode("indexes")

extract_col = ['peptide_sequence','modeSequence','phosphosite','cleaned_site',
                'locProbability','p_value','avg_ratio_absfc','log2fc_log10fc','exp_condition_old',
                'mapped_phosphosite',
                'expression','exp_cond_id']

newdf = pd.DataFrame(columns=df.columns)

for _, row in df.iterrows():
	li_in = row['indexes']
	row["peptide_sequence"] = [row["peptide_sequence"][i] for i in li_in]
	row["modeSequence"] = [row["modeSequence"][i] for i in li_in]
	row["phosphosite"] = [row["phosphosite"][i] for i in li_in]
	row["cleaned_site"] = [row["cleaned_site"][i] for i in li_in]
	row["locProbability"] = [row["locProbability"][i] for i in li_in]
	row["p_value"] = [row["p_value"][i] for i in li_in]
	row["avg_ratio_absfc"] = [row["avg_ratio_absfc"][i] for i in li_in]
	row["log2fc_log10fc"] = [row["log2fc_log10fc"][i] for i in li_in]
	row["exp_condition_old"] = [row["exp_condition_old"][i] for i in li_in]
	row["mapped_phosphosite"] = [row["mapped_phosphosite"][i] for i in li_in]
	row["expression"] = [row["expression"][i] for i in li_in]
	row["exp_cond_id"] = [row["exp_cond_id"][i] for i in li_in]
	newdf = newdf.append(row, ignore_index=True)

newdf['checksame'] = newdf['expression'].apply(lambda x:True if len(set(x)) == 1 else False)

same_df = newdf.loc[ newdf["checksame"] == True ]
same_df = same_df[final_cols]

FINAL_dF = pd.concat([FINAL_dF,same_df])

df = newdf.loc[ newdf["checksame"] == False ]

newdf_2 = pd.DataFrame()
for index_, row in df.iterrows():
	if row["cleaned_site"].count(row["mapped_phosphosite"][0]) == 1:
		i = row["cleaned_site"].index(row["mapped_phosphosite"][0])
		print(row["expression"][i])
		row["peptide_sequence"] = [row["peptide_sequence"][i]]
		row["modeSequence"] = [row["modeSequence"][i]]
		row["phosphosite"] = [row["phosphosite"][i] ]
		row["cleaned_site"] = [row["cleaned_site"][i] ]
		row["locProbability"] = [row["locProbability"][i]]
		row["p_value"] = [row["p_value"][i]]
		row["avg_ratio_absfc"] = [row["avg_ratio_absfc"][i] ]
		row["log2fc_log10fc"] = [row["log2fc_log10fc"][i] ]
		row["exp_condition_old"] = [row["exp_condition_old"][i] ]
		row["mapped_phosphosite"] = [row["mapped_phosphosite"][i] ]
		row["expression"] = [row["expression"][i] ]
		row["exp_cond_id"] = [row["exp_cond_id"][i] ]
		newdf_2 = newdf_2.append(row, ignore_index=True)
	else:
		print(index_)
		newdf_2 = newdf_2.append(row, ignore_index=True)

newdf_2["mappedd"] = newdf_2["expression"].apply(lambda x:True if len(set(x)) == 1 else False )

mapped_df = newdf_2.loc[newdf_2["mappedd"] == True ]

mapped_df = mapped_df[final_cols]
FINAL_dF = pd.concat([FINAL_dF,mapped_df])


unmapped_df = newdf_2.loc[newdf_2["mappedd"] == False ]
unmapped_df.drop("mappedd", axis=1, inplace=True)

def getmax_expression(x):
	up_c = x.count("Up-regulated")
	down_c = x.count("down-regulated")

	if up_c == down_c:
		return x
	elif up_c > down_c:
		return ["Up-regulated"]
	else:
		return ["down-regulated"]

unmapped_df["expression"] = unmapped_df["expression"].apply(getmax_expression)
unmapped_df["mappedd"] = unmapped_df["expression"].apply(lambda x:True if len(set(x)) == 1 else False )

mapped_df = unmapped_df.loc[unmapped_df["mappedd"] == True ]

FInal_unmapped = unmapped_df.loc[unmapped_df["mappedd"] == False ]
FInal_unmapped.to_excel('FInal_unmapped.xlsx')

mapped_df = mapped_df[final_cols]

FINAL_dF = pd.concat([FINAL_dF,mapped_df])

FINAL_dF["mapped_phosphosite"] = FINAL_dF["mapped_phosphosite"].apply(lambda x:list(set(x)))
FINAL_dF["expression"] = FINAL_dF["expression"].apply(lambda x:list(set(x)))

FINAL_dF.to_excel('FInal_MAPPED.xlsx')








