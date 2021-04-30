import pandas as pd
import numpy as np
import collections
import matplotlib.pyplot as plt
import math

failinimede_list = ["Alasoo_2018.macrophage_IFNg+Salmonella_", "Alasoo_2018.macrophage_IFNg_",
                    "Alasoo_2018.macrophage_naive_", 'BLUEPRINT_PE.T-cell_', 'BLUEPRINT_SE.monocyte_',
                    'BLUEPRINT_SE.neutrophil_', 'BrainSeq.brain_', 'FUSION.adipose_naive_', 'FUSION.muscle_naive_',
                    'GENCORD.fibroblast_', 'GENCORD.LCL_', 'GENCORD.T-cell_', 'GEUVADIS.LCL_',
                    'GTEx.adipose_subcutaneous_', 'GTEx.adipose_visceral_', 'GTEx.adrenal_gland_', 'GTEx.artery_aorta_',
                    'GTEx.artery_coronary_', 'GTEx.artery_tibial_', 'GTEx.blood_', 'GTEx.brain_amygdala_',
                    'GTEx.brain_anterior_cingulate_cortex_', 'GTEx.brain_caudate_', 'GTEx.brain_cerebellar_hemisphere_',
                    'GTEx.brain_cerebellum_', 'GTEx.brain_cortex_', 'GTEx.brain_frontal_cortex_',
                    'GTEx.brain_hippocampus_', 'GTEx.brain_hypothalamus_', 'GTEx.brain_nucleus_accumbens_',
                    'GTEx.brain_putamen_', 'GTEx.brain_spinal_cord_', 'GTEx.brain_substantia_nigra_', 'GTEx.breast_',
                    'GTEx.colon_sigmoid_', 'GTEx.colon_transverse_', 'GTEx.esophagus_gej_', 'GTEx.esophagus_mucosa_',
                    'GTEx.esophagus_muscularis_', 'GTEx.fibroblast_', 'GTEx.heart_atrial_appendage_',
                    'GTEx.heart_left_ventricle_', 'GTEx.LCL_', 'GTEx.liver_', 'GTEx.lung_',
                    'GTEx.minor_salivary_gland_', 'GTEx.muscle_', 'GTEx.nerve_tibial_', 'GTEx.ovary_', 'GTEx.pancreas_',
                    'GTEx.pituitary_', 'GTEx.prostate_', 'GTEx.skin_not_sun_exposed_', 'GTEx.skin_sun_exposed_',
                    'GTEx.small_intestine_', 'GTEx.spleen_', 'GTEx.stomach_', 'GTEx.testis_', 'GTEx.thyroid_',
                    'GTEx.uterus_', 'GTEx.vagina_', 'HipSci.iPSC_', 'Lepik_2017.blood_',
                    'Nedelec_2016.macrophage_Listeria_', 'Nedelec_2016.macrophage_naive_',
                    'Nedelec_2016.macrophage_Salmonella_', 'Quach_2016.monocyte_IAV_', 'Quach_2016.monocyte_LPS_',
                    'Quach_2016.monocyte_naive_', 'Quach_2016.monocyte_Pam3CSK4_', 'Quach_2016.monocyte_R848_',
                    'ROSMAP.brain_naive_', 'Schmiedel_2018.B-cell_naive_', 'Schmiedel_2018.CD4_T-cell_anti-CD3-CD28_',
                    'Schmiedel_2018.CD4_T-cell_naive_', 'Schmiedel_2018.CD8_T-cell_anti-CD3-CD28_',
                    'Schmiedel_2018.CD8_T-cell_naive_', 'Schmiedel_2018.monocyte_CD16_naive_',
                    'Schmiedel_2018.monocyte_naive_', 'Schmiedel_2018.NK-cell_naive_', 'Schmiedel_2018.Tfh_memory_',
                    'Schmiedel_2018.Th1-17_memory_', 'Schmiedel_2018.Th17_memory_', 'Schmiedel_2018.Th1_memory_',
                    'Schmiedel_2018.Th2_memory_', 'Schmiedel_2018.Treg_memory_', 'Schmiedel_2018.Treg_naive_',
                    'Schwartzentruber_2018.sensory_neuron_', 'TwinsUK.blood_', 'TwinsUK.fat_', 'TwinsUK.LCL_',
                    'TwinsUK.skin_', 'van_de_Bunt_2015.pancreatic_islet_']


# failinimede_list = ["HipSci.iPSC_", "Alasoo_2018.macrophage_IFNg+Salmonella_", "Alasoo_2018.macrophage_IFNg_", "Alasoo_2018.macrophage_naive_"]


def ümarda(number):
    return int(math.ceil(number / 10.0)) * 10


def leia_cs_suurused(failinimede_list, loendus_kuni, tüüp):
    cs_suurused = {x: 0 for x in range(10, loendus_kuni + 11, 10)}

    if tüüp == 1:  # geeniekspressioon
        laiend = "ge"
    elif tüüp == 2:  # RNA splaissimine
        laiend = "txrev"
    elif tüüp == 3:  # promootor
        laiend = "txrev"
    elif tüüp == 4:
        laiend = "txrev"
    elif tüüp == 5:
        laiend = "tx"
    else:
        laiend = "muu"

    lisa = laiend + ".purity_filtered.txt"

    for fail in failinimede_list:
        tabel = pd.read_csv(r"C:\Users\pr\Documents\kool\6 semester 2020-21\Baka\andmed\eraldi\\" + fail + lisa,
                            sep="	")  # loeb faili sisse
        tabel = tabel[['phenotype_id', 'cs_id', 'cs_size', 'variant_id']]

        if tüüp == 2:
            tabel = tabel[tabel['phenotype_id'].str.contains("contained")]  # jätab alles vaid andmed splaissimise kohta
        elif tüüp == 3:
            tabel = tabel[tabel['phenotype_id'].str.contains("upstream")]
        elif tüüp == 4:
            tabel = tabel[tabel['phenotype_id'].str.contains("downstream")]

        tabel = tabel.drop_duplicates(["cs_id"])  # jätame cs-id väid ühekordselt alles
        for i, row in tabel.iterrows():
            if (row.cs_size > loendus_kuni):  # ületab etteantud piiri (100)
                cs_suurused[loendus_kuni + 10] += 1
            else:
                cs_suurused[ümarda(row.cs_size)] += 1

    cs_suurused = collections.OrderedDict(sorted(cs_suurused.items()))  # järjestab sõnastiku
    print(cs_suurused)

    return cs_suurused


loendus_kuni = 100
ge_cs_suurused = leia_cs_suurused(failinimede_list, loendus_kuni, 1)
txrev_cs_suurused = leia_cs_suurused(failinimede_list, loendus_kuni, 2)
upstream_cs_suurused = leia_cs_suurused(failinimede_list, loendus_kuni, 3)
downstream_cs_suurused = leia_cs_suurused(failinimede_list, loendus_kuni, 4)
tx_cs_suurused = leia_cs_suurused(failinimede_list, loendus_kuni, 5)

print("Geeniekspressioom")
print(ge_cs_suurused)

print("\nRNA spl")
print(txrev_cs_suurused)

print("\nupstream")
print(upstream_cs_suurused)

print("\ndownstream")
print(downstream_cs_suurused)

print("\ntranscript")
print(tx_cs_suurused)

labels = ["1-10", "11-20", "21-30", "31-40", "41-50", "51-60", "61-70", "71-80", "81-90", "91-100", "101+"]

x = np.arange(len(labels))  # the label locations
bar_width = 0.15  # the width of the bars

fig, ax = plt.subplots(figsize=(8, 4))

bars1 = list(ge_cs_suurused.values())
bars2 = list(txrev_cs_suurused.values())
bars3 = list(upstream_cs_suurused.values())
bars4 = list(downstream_cs_suurused.values())
bars5 = list(tx_cs_suurused.values())

# Set position of bar on X axis
r1 = np.arange(len(bars1))
r2 = [x + bar_width for x in r1]
r3 = [x + bar_width for x in r2]
r4 = [x + bar_width for x in r3]
r5 = [x + bar_width for x in r4]

# Make the plot
ax.bar(r1, bars1, color='c', width=bar_width, edgecolor='white', label='Geeniekspressioon')
ax.bar(r2, bars2, color='salmon', width=bar_width, edgecolor='white', label='RNA splaissimine')
ax.bar(r3, bars3, color='orchid', width=bar_width, edgecolor='white', label='Promootor')
ax.bar(r4, bars4, color='lightgreen', width=bar_width, edgecolor='white', label='Terminaator')
ax.bar(r5, bars5, color='gold', width=bar_width, edgecolor='white', label='Transkript')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Leitud CS-ide arv')
ax.set_xlabel('CS suurus')
ax.set_title('CS-ide suurused vastavalt uurimistüübile')
ax.set_xticks(x)

ax.set_xticklabels(labels)
ax.legend()

fig.tight_layout()

plt.show()
