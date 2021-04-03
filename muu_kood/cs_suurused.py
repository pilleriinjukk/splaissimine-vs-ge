import pandas as pd
import numpy as np
from collections import defaultdict
import collections
from matplotlib.axes import Axes
import matplotlib.pyplot as plt
import math



ge_failinimede_list = ["Alasoo_2018.macrophage_IFNg+Salmonella_ge.purity_filtered.txt", "Alasoo_2018.macrophage_IFNg_ge.purity_filtered.txt", "Alasoo_2018.macrophage_naive_ge.purity_filtered.txt", "Alasoo_2018.macrophage_Salmonella_ge.purity_filtered.txt", "BLUEPRINT_PE.T-cell_ge.purity_filtered.txt", "BLUEPRINT_SE.monocyte_ge.purity_filtered.txt", "BLUEPRINT_SE.neutrophil_ge.purity_filtered.txt", "BrainSeq.brain_ge.purity_filtered.txt", "FUSION.adipose_naive_ge.purity_filtered.txt", "FUSION.muscle_naive_ge.purity_filtered.txt", "GENCORD.fibroblast_ge.purity_filtered.txt", "GENCORD.LCL_ge.purity_filtered.txt", "GENCORD.T-cell_ge.purity_filtered.txt", "GEUVADIS.LCL_ge.purity_filtered.txt", "GTEx.adipose_subcutaneous_ge.purity_filtered.txt", "GTEx.adipose_visceral_ge.purity_filtered.txt", "GTEx.adrenal_gland_ge.purity_filtered.txt", "GTEx.artery_aorta_ge.purity_filtered.txt", "GTEx.artery_coronary_ge.purity_filtered.txt", "GTEx.artery_tibial_ge.purity_filtered.txt", "GTEx.blood_ge.purity_filtered.txt", "GTEx.brain_amygdala_ge.purity_filtered.txt", "GTEx.brain_anterior_cingulate_cortex_ge.purity_filtered.txt", "GTEx.brain_caudate_ge.purity_filtered.txt", "GTEx.brain_cerebellar_hemisphere_ge.purity_filtered.txt", "GTEx.brain_cerebellum_ge.purity_filtered.txt", "GTEx.brain_cortex_ge.purity_filtered.txt", "GTEx.brain_frontal_cortex_ge.purity_filtered.txt", "GTEx.brain_hippocampus_ge.purity_filtered.txt", "GTEx.brain_hypothalamus_ge.purity_filtered.txt", "GTEx.brain_nucleus_accumbens_ge.purity_filtered.txt", "GTEx.brain_putamen_ge.purity_filtered.txt", "GTEx.brain_spinal_cord_ge.purity_filtered.txt", "GTEx.brain_substantia_nigra_ge.purity_filtered.txt", "GTEx.breast_ge.purity_filtered.txt", "GTEx.colon_sigmoid_ge.purity_filtered.txt", "GTEx.colon_transverse_ge.purity_filtered.txt", "GTEx.esophagus_gej_ge.purity_filtered.txt", "GTEx.esophagus_mucosa_ge.purity_filtered.txt", "GTEx.esophagus_muscularis_ge.purity_filtered.txt", "GTEx.fibroblast_ge.purity_filtered.txt", "GTEx.heart_atrial_appendage_ge.purity_filtered.txt", "GTEx.heart_left_ventricle_ge.purity_filtered.txt", "GTEx.LCL_ge.purity_filtered.txt", "GTEx.liver_ge.purity_filtered.txt", "GTEx.lung_ge.purity_filtered.txt", "GTEx.minor_salivary_gland_ge.purity_filtered.txt", "GTEx.muscle_ge.purity_filtered.txt", "GTEx.nerve_tibial_ge.purity_filtered.txt", "GTEx.ovary_ge.purity_filtered.txt", "GTEx.pancreas_ge.purity_filtered.txt", "GTEx.pituitary_ge.purity_filtered.txt", "GTEx.prostate_ge.purity_filtered.txt", "GTEx.skin_not_sun_exposed_ge.purity_filtered.txt", "GTEx.skin_sun_exposed_ge.purity_filtered.txt", "GTEx.small_intestine_ge.purity_filtered.txt", "GTEx.spleen_ge.purity_filtered.txt", "GTEx.stomach_ge.purity_filtered.txt", "GTEx.testis_ge.purity_filtered.txt", "GTEx.thyroid_ge.purity_filtered.txt", "GTEx.uterus_ge.purity_filtered.txt", "GTEx.vagina_ge.purity_filtered.txt", "HipSci.iPSC_ge.purity_filtered.txt", "Lepik_2017.blood_ge.purity_filtered.txt", "Nedelec_2016.macrophage_Listeria_ge.purity_filtered.txt", "Nedelec_2016.macrophage_naive_ge.purity_filtered.txt", "Nedelec_2016.macrophage_Salmonella_ge.purity_filtered.txt", "Quach_2016.monocyte_IAV_ge.purity_filtered.txt", "Quach_2016.monocyte_LPS_ge.purity_filtered.txt", "Quach_2016.monocyte_naive_ge.purity_filtered.txt", "Quach_2016.monocyte_Pam3CSK4_ge.purity_filtered.txt", "Quach_2016.monocyte_R848_ge.purity_filtered.txt", "ROSMAP.brain_naive_ge.purity_filtered.txt", "Schmiedel_2018.B-cell_naive_ge.purity_filtered.txt", "Schmiedel_2018.CD4_T-cell_anti-CD3-CD28_ge.purity_filtered.txt", "Schmiedel_2018.CD4_T-cell_naive_ge.purity_filtered.txt", "Schmiedel_2018.CD8_T-cell_anti-CD3-CD28_ge.purity_filtered.txt", "Schmiedel_2018.CD8_T-cell_naive_ge.purity_filtered.txt", "Schmiedel_2018.monocyte_CD16_naive_ge.purity_filtered.txt", "Schmiedel_2018.monocyte_naive_ge.purity_filtered.txt", "Schmiedel_2018.NK-cell_naive_ge.purity_filtered.txt", "Schmiedel_2018.Tfh_memory_ge.purity_filtered.txt", "Schmiedel_2018.Th1-17_memory_ge.purity_filtered.txt", "Schmiedel_2018.Th17_memory_ge.purity_filtered.txt", "Schmiedel_2018.Th1_memory_ge.purity_filtered.txt", "Schmiedel_2018.Th2_memory_ge.purity_filtered.txt", "Schmiedel_2018.Treg_memory_ge.purity_filtered.txt", "Schmiedel_2018.Treg_naive_ge.purity_filtered.txt", "Schwartzentruber_2018.sensory_neuron_ge.purity_filtered.txt", "TwinsUK.blood_ge.purity_filtered.txt", "TwinsUK.fat_ge.purity_filtered.txt", "TwinsUK.LCL_ge.purity_filtered.txt", "TwinsUK.skin_ge.purity_filtered.txt", "van_de_Bunt_2015.pancreatic_islet_ge.purity_filtered.txt"]
txrev_failinimede_list = ["Alasoo_2018.macrophage_IFNg+Salmonella_txrev.purity_filtered.txt", "Alasoo_2018.macrophage_IFNg_txrev.purity_filtered.txt", "Alasoo_2018.macrophage_naive_txrev.purity_filtered.txt", "Alasoo_2018.macrophage_Salmonella_txrev.purity_filtered.txt", "BLUEPRINT_PE.T-cell_txrev.purity_filtered.txt", "BLUEPRINT_SE.monocyte_txrev.purity_filtered.txt", "BLUEPRINT_SE.neutrophil_txrev.purity_filtered.txt", "BrainSeq.brain_txrev.purity_filtered.txt", "FUSION.adipose_naive_txrev.purity_filtered.txt", "FUSION.muscle_naive_txrev.purity_filtered.txt", "GENCORD.fibroblast_txrev.purity_filtered.txt", "GENCORD.LCL_txrev.purity_filtered.txt", "GENCORD.T-cell_txrev.purity_filtered.txt", "GEUVADIS.LCL_txrev.purity_filtered.txt", "GTEx.adipose_subcutaneous_txrev.purity_filtered.txt", "GTEx.adipose_visceral_txrev.purity_filtered.txt", "GTEx.adrenal_gland_txrev.purity_filtered.txt", "GTEx.artery_aorta_txrev.purity_filtered.txt", "GTEx.artery_coronary_txrev.purity_filtered.txt", "GTEx.artery_tibial_txrev.purity_filtered.txt", "GTEx.blood_txrev.purity_filtered.txt", "GTEx.brain_amygdala_txrev.purity_filtered.txt", "GTEx.brain_anterior_cingulate_cortex_txrev.purity_filtered.txt", "GTEx.brain_caudate_txrev.purity_filtered.txt", "GTEx.brain_cerebellar_hemisphere_txrev.purity_filtered.txt", "GTEx.brain_cerebellum_txrev.purity_filtered.txt", "GTEx.brain_cortex_txrev.purity_filtered.txt", "GTEx.brain_frontal_cortex_txrev.purity_filtered.txt", "GTEx.brain_hippocampus_txrev.purity_filtered.txt", "GTEx.brain_hypothalamus_txrev.purity_filtered.txt", "GTEx.brain_nucleus_accumbens_txrev.purity_filtered.txt", "GTEx.brain_putamen_txrev.purity_filtered.txt", "GTEx.brain_spinal_cord_txrev.purity_filtered.txt", "GTEx.brain_substantia_nigra_txrev.purity_filtered.txt", "GTEx.breast_txrev.purity_filtered.txt", "GTEx.colon_sigmoid_txrev.purity_filtered.txt", "GTEx.colon_transverse_txrev.purity_filtered.txt", "GTEx.esophagus_gej_txrev.purity_filtered.txt", "GTEx.esophagus_mucosa_txrev.purity_filtered.txt", "GTEx.esophagus_muscularis_txrev.purity_filtered.txt", "GTEx.fibroblast_txrev.purity_filtered.txt", "GTEx.heart_atrial_appendage_txrev.purity_filtered.txt", "GTEx.heart_left_ventricle_txrev.purity_filtered.txt", "GTEx.LCL_txrev.purity_filtered.txt", "GTEx.liver_txrev.purity_filtered.txt", "GTEx.lung_txrev.purity_filtered.txt", "GTEx.minor_salivary_gland_txrev.purity_filtered.txt", "GTEx.muscle_txrev.purity_filtered.txt", "GTEx.nerve_tibial_txrev.purity_filtered.txt", "GTEx.ovary_txrev.purity_filtered.txt", "GTEx.pancreas_txrev.purity_filtered.txt", "GTEx.pituitary_txrev.purity_filtered.txt", "GTEx.prostate_txrev.purity_filtered.txt", "GTEx.skin_not_sun_exposed_txrev.purity_filtered.txt", "GTEx.skin_sun_exposed_txrev.purity_filtered.txt", "GTEx.small_intestine_txrev.purity_filtered.txt", "GTEx.spleen_txrev.purity_filtered.txt", "GTEx.stomach_txrev.purity_filtered.txt", "GTEx.testis_txrev.purity_filtered.txt", "GTEx.thyroid_txrev.purity_filtered.txt", "GTEx.uterus_txrev.purity_filtered.txt", "GTEx.vagina_txrev.purity_filtered.txt", "HipSci.iPSC_txrev.purity_filtered.txt", "Lepik_2017.blood_txrev.purity_filtered.txt", "Nedelec_2016.macrophage_Listeria_txrev.purity_filtered.txt", "Nedelec_2016.macrophage_naive_txrev.purity_filtered.txt", "Nedelec_2016.macrophage_Salmonella_txrev.purity_filtered.txt", "Quach_2016.monocyte_IAV_txrev.purity_filtered.txt", "Quach_2016.monocyte_LPS_txrev.purity_filtered.txt", "Quach_2016.monocyte_naive_txrev.purity_filtered.txt", "Quach_2016.monocyte_Pam3CSK4_txrev.purity_filtered.txt", "Quach_2016.monocyte_R848_txrev.purity_filtered.txt", "ROSMAP.brain_naive_txrev.purity_filtered.txt", "Schmiedel_2018.B-cell_naive_txrev.purity_filtered.txt", "Schmiedel_2018.CD4_T-cell_anti-CD3-CD28_txrev.purity_filtered.txt", "Schmiedel_2018.CD4_T-cell_naive_txrev.purity_filtered.txt", "Schmiedel_2018.CD8_T-cell_anti-CD3-CD28_txrev.purity_filtered.txt", "Schmiedel_2018.CD8_T-cell_naive_txrev.purity_filtered.txt", "Schmiedel_2018.monocyte_CD16_naive_txrev.purity_filtered.txt", "Schmiedel_2018.monocyte_naive_txrev.purity_filtered.txt", "Schmiedel_2018.NK-cell_naive_txrev.purity_filtered.txt", "Schmiedel_2018.Tfh_memory_txrev.purity_filtered.txt", "Schmiedel_2018.Th1-17_memory_txrev.purity_filtered.txt", "Schmiedel_2018.Th17_memory_txrev.purity_filtered.txt", "Schmiedel_2018.Th1_memory_txrev.purity_filtered.txt", "Schmiedel_2018.Th2_memory_txrev.purity_filtered.txt", "Schmiedel_2018.Treg_memory_txrev.purity_filtered.txt", "Schmiedel_2018.Treg_naive_txrev.purity_filtered.txt", "Schwartzentruber_2018.sensory_neuron_txrev.purity_filtered.txt", "TwinsUK.blood_txrev.purity_filtered.txt", "TwinsUK.fat_txrev.purity_filtered.txt", "TwinsUK.LCL_txrev.purity_filtered.txt", "TwinsUK.skin_txrev.purity_filtered.txt", "van_de_Bunt_2015.pancreatic_islet_txrev.purity_filtered.txt"]
# ge_failinimede_list = ["GTEx.minor_salivary_gland_ge.purity_filtered.txt"]
# txrev_failinimede_list = ["GTEx.minor_salivary_gland_txrev.purity_filtered.txt"]




def ümarda(number):
    return int(math.ceil(number / 10.0)) * 10

def leia_cs_suurused(failinimede_list, loendus_kuni, tüüp):
    cs_suurused = {x: 0 for x in range(10, loendus_kuni + 1, 10)}

    for fail in failinimede_list:
        tabel = pd.read_csv(r"C:\Users\pr\Documents\kool\6 semester 2020-21\Baka\andmed\\" + fail, sep="	")  # loeb faili sisse
        tabel = tabel[['phenotype_id', 'cs_id', 'cs_size', 'variant_id']]

        #kui on splaissimine
        if tüüp == 2:
            tabel = tabel[tabel['phenotype_id'].str.contains("contained")]  # jätab alles vaid andmed splaissimise kohta
        elif tüüp == 3:
            tabel = tabel[tabel['phenotype_id'].str.contains("upstream")]
        elif tüüp == 4:
            tabel = tabel[tabel['phenotype_id'].str.contains("downstream")]

        tabel = tabel.drop_duplicates(["cs_id"])
        for i, row in tabel.iterrows():
            if (row.cs_size > loendus_kuni):
                cs_suurused[loendus_kuni] += 1
            else:
                cs_suurused[ümarda(row.cs_size)] += 1

    #cs_suurused = collections.OrderedDict(sorted(cs_suurused.items()))  # järjestab sõnastiku
    print(cs_suurused)

    return cs_suurused

loendus_kuni = 100
ge_cs_suurused = leia_cs_suurused(ge_failinimede_list, loendus_kuni, 1)
txrev_cs_suurused = leia_cs_suurused(txrev_failinimede_list, loendus_kuni, 2)
upstream_cs_suurused = leia_cs_suurused(txrev_failinimede_list, loendus_kuni, 3)
downstream_cs_suurused = leia_cs_suurused(txrev_failinimede_list, loendus_kuni, 4)



labels = list(range(10, loendus_kuni + 1, 10))


#plt.bar(cs_suurused.keys(), cs_suurused.values(), color="lightgreen")

x = np.arange(len(labels))  # the label locations
barWidth = 0.2  # the width of the bars

fig, ax = plt.subplots(figsize=(8, 4))

bars1 = list(ge_cs_suurused.values())
bars2 = list(txrev_cs_suurused.values())
bars3 = list(upstream_cs_suurused.values())
bars4 = list(downstream_cs_suurused.values())


# Set position of bar on X axis
r1 = np.arange(len(bars1))
r2 = [x + barWidth for x in r1]
r3 = [x + barWidth for x in r2]
r4 = [x + barWidth for x in r3]

# Make the plot
ax.bar(r1, bars1, color='c', width=barWidth, edgecolor='white', label='Geeniekspressioon')
ax.bar(r2, bars2, color='salmon', width=barWidth, edgecolor='white', label='RNA splaissimine')
ax.bar(r3, bars3, color='orchid', width=barWidth, edgecolor='white', label='Upstream')
ax.bar(r4, bars4, color='lightgreen', width=barWidth, edgecolor='white', label='Downstream')


# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Leitud CS-ide arv')
ax.set_xlabel('CS suurus')
ax.set_title('Põhjuslike variantide kogumite (CS-ide) suurused vastavalt uurimistüübile')
ax.set_xticks(x)
labels = labels[:-1]
labels.append("+" + str(loendus_kuni))
ax.set_xticklabels(labels)
ax.legend()


#ax.Axes.bar_label(rects1, padding=3)
#ax.Axes.bar_label(rects2, padding=3)

fig.tight_layout()

plt.show()

