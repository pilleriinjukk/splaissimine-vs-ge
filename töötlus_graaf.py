import pandas as pd
import pyranges as pr
import numpy as np
import matplotlib.pyplot as plt
import collections
from collections import defaultdict
from os import path
import igraph as ig
from igraph import VertexSeq
from itertools import combinations


pd.set_option('display.max_columns', None)
desired_width = 1000
pd.set_option('display.width', desired_width)
pd.set_option('display.max_colwidth', 100)
np.set_printoptions(linewidth=desired_width)
pd.options.mode.chained_assignment = None

'''
Failist andmete sisselugemine, veergude lisamine ja kustutamine ning kromosoomipõhiselt failidesse salvestamine.
Sisend: failinimi, sisseloetavate andmete tüüp (geeniekspressioon, splaissimine, promootor)
'''
def algtöötlus(failinimi, tüüp=""):
    tabel = pd.read_csv(r"C:\Users\pr\Documents\kool\6 semester 2020-21\Baka\andmed\\" + failinimi, sep="	")  # loeb faili sisse
    tabel = tabel[['phenotype_id', 'chr', 'cs_id', 'pip', 'cs_size', 'variant_id']]  # mis veerud jätta alles
    tabel = tabel.rename(columns={'chr': 'Chromosome'})
    tabel = pd.concat([tabel, tabel["variant_id"].str.split("_", expand=True)[1].rename("Start")], axis=1)  # uus veerg Start
    tabel["Start"] = pd.to_numeric(tabel["Start"])
    tabel['End'] = tabel['Start'] + 1  # Pyranges puhul peab mingi vahemik olema, suurema vahemiku puhul loetakse liiga palju ülekatteid
    tabel = tabel[['phenotype_id', 'cs_id', 'Chromosome', 'Start', 'End', 'pip', 'cs_size']]  # mis veerud jätta alles

    nimi = "ge"
    if tüüp == "spalissimine":  # vaid txrev puhul
        tabel = tabel[tabel['phenotype_id'].str.contains("contained")]  # jätab alles vaid andmed splaissimise kohta
        nimi = "txrev"

    if tüüp == "promootor":  # vaid txrev puhul
        tabel = tabel[tabel['phenotype_id'].str.contains("upstream")]  # jätab alles vaid andmed splaissimise kohta
        nimi = "txrev_up"

    tabel["cs_id"] = tabel["cs_id"] + "_" + failinimi.split(".purity_filtered.txt")[0] #lisame cs_id veerule andmebaasi nime juurde

    #teeme igale kromosoomile oma tabeli
    tabelid = list(range(1)) #TODO muidu 22

    for i in range(1): #kirjutame iga kromosoomi oma faili #TODO muidu 22
        tabelid[i] = (tabel[(tabel['Chromosome'] == 1 + i)])
        if path.exists(r"C:\Users\pr\Documents\kool\6 semester 2020-21\Baka\andmed\\" + nimi + "_" + str(i + 1) + ".csv"):
            tabelid[i].to_csv(
                r"C:\Users\pr\Documents\kool\6 semester 2020-21\Baka\andmed\\" + nimi + "_" + str(i + 1) + ".csv", mode='a', index=False, header=False)
        else:
            tabelid[i].to_csv(
                r"C:\Users\pr\Documents\kool\6 semester 2020-21\Baka\andmed\\" + nimi + "_" + str(i + 1) + ".csv", mode='a', index=False)  # esimene kord salvestame koos veerunimetustega


'''
Tabelist ülekatete leidmine kasutades Pyranges klasterdamist. Joonistab klastrid, kus tippudeks on cs-id ja servad on nende vahel siis kui need jagavad omavahel samu variante
Sisend: tabel, faili nime algus, maksimaalne credible set-i suurus, kas on tegemist splaissimisega
Väljund: leitud ülekatted
'''
def ülekatete_leidmine(tabel, nimi, cs_suurus, on_splaissimine):
    tabel = tabel[tabel['cs_size'] <= cs_suurus]  # credible set maksimaalne suurus
    if on_splaissimine:
        tabel['phenotype_id'] = tabel['phenotype_id'].str.split(".", expand=True)  # jätab alles vaid geeni nime

    pr_tabel = pr.PyRanges(tabel)  # teeb Pyranges tabeli
    cluster = pr_tabel.cluster(by="Start") #klasterdame asukohapõhiselt
    cluster_df = cluster.df
    # cluster_df = cluster_df[:50] #TODO: eemalda see

    # jätame alles vaid sellised klastrid, kus on vähemalt kaks cs-i
    cluster_values = cluster_df.Cluster.value_counts()  # type - series
    cluster_values = cluster_values[cluster_values > 1]
    cluster_id = cluster_values.axes[0].tolist()
    for i, rida in cluster_df.iterrows():
        if (rida["Cluster"] not in cluster_id):
            cluster_df = cluster_df.drop(index=i)

    #saame kätte cs-id
    cs_id = set((cluster_df.cs_id.values))

    g = ig.Graph()

    # tippude lisamine
    g.add_vertices(len(cs_id) + 1)
    g.vs[-1]['cs'] = "lisa"
    g.vs[-1]["cluster"] = []
    g.vs[-1]["phenotype_id"] = "lisa"
    g.vs[-1]["positions"] = []

    lisatava_tipu_indeks = 0
    eelmine_cluster = -1
    lisa_servad = []
    for i, rida in cluster_df.iterrows():
        cluster = rida.Cluster
        #lisame tipu
        try: #kui leidub tipp
            indeks = g.vs.find(cs=rida.cs_id).index
            g.vs[indeks]['cluster'].append(cluster)
            g.vs[indeks]["positions"].append(rida.Start)
        except: #sellist tippu veel pole
            indeks = lisatava_tipu_indeks
            lisatava_tipu_indeks += 1
            g.vs[indeks]['cs'] = rida.cs_id
            g.vs[indeks]['cluster'] = [cluster]
            g.vs[indeks]["positions"] = [rida.Start]
            g.vs[indeks]["phenotype_id"] = rida.phenotype_id

        #servad
        if (cluster == eelmine_cluster):
            lisa_servad.append(indeks)
        else:
            #teeme servad eelmiste klastritega
            for servad in combinations(lisa_servad, 2):
                g.add_edge(servad[0], servad[1])

            #muudame praeguse klastri informatsiooni
            g.summary()
            lisa_servad = [indeks]
            eelmine_cluster = cluster

    for servad in combinations(lisa_servad, 2):
        g.add_edge(servad[0], servad[1])

    g.delete_vertices(len(cs_id)) #eemadlame lisa rea


    kustutamisele = []
    components = g.clusters()
    geenid_ülekattes = []
    komponentide_suurused = defaultdict(int)
    komponentide_geenid = []

    #leiame, mis geenid kuuluvad samadesse usalduskomponenti
    for component in components:
        geenid = set()
        for vertex in component:
            geenid.add(g.vs[vertex]["phenotype_id"])

        # if (len(geenid) == 1): #cs-id on ülekattes vaid üle ühe geeni
        #     for vertex in component:
        #         kustutamisele.append(vertex)
        # else:

        komponentide_suurused[len(component)] += 1
        komponentide_geenid.append(geenid)

        geenid_ülekattes.append(geenid)

    print("Komponentide suurused:", komponentide_suurused)
    print("Ühte komponenti kuuluvad geenid:", komponentide_geenid)

    #eemaldame klastrid, kus olid ülekatted vaid üle ühe geeni
    #g.delete_vertices(kustutamisele)

    #salvestame leitud ülekatted
    ülekatteid = defaultdict(int)
    for komponent in geenid_ülekattes:
        ülekatteid[len(komponent)] += 1

    #salvestame tipud faili
    tipud_df = g.get_vertex_dataframe()
    #print(tipud_df)
    tipud_df.to_csv(nimi + "_tipud.csv")

    #joonistame
    components = g.clusters()
    layout = g.layout(layout='auto')
    ig.plot(components)

    return ülekatteid


'''
Ülekatete sõnastike sama pikkuseks tegemine ja järjestamine, et oleks lihtsam joonist teha
Sisend: ge ja txrev ülekatete sõnastikud
Väljund: korrastatud ge ja txrev sõnastikud
'''
def ülekatete_korrastamine(ge_count, txrev_count, ge_tabelite_pikkus, txrev_tabelite_pikkus):
    suurim_pikkus = max(len(ge_count), len(txrev_count))

    for sõnastik in (ge_count, txrev_count):
        while (len(sõnastik) < suurim_pikkus):
            sõnastik[len(sõnastik) + 1] = 0

    ge_count = collections.OrderedDict(sorted(ge_count.items()))  # järjestab sõnastiku
    txrev_count = collections.OrderedDict(sorted(txrev_count.items()))  # järjestab sõnastiku

    ge_labels = list(ge_count.values())
    txrev_labels = list(txrev_count.values())

    for count in ge_count:
        ge_count[count] = ge_count[count] / ge_tabelite_pikkus

    for count in txrev_count:
        txrev_count[count] = txrev_count[count] / txrev_tabelite_pikkus

    return ge_count, txrev_count, ge_labels, txrev_labels


cs_suurus = 10  # credible set maksimaalne suurus

# geeniekspressioon
print("\nGeeniekspressioon:")
ge_failinimede_list = ["Alasoo_2018.macrophage_IFNg+Salmonella_ge.purity_filtered.txt", "Alasoo_2018.macrophage_IFNg_ge.purity_filtered.txt", "Alasoo_2018.macrophage_naive_ge.purity_filtered.txt", "Alasoo_2018.macrophage_Salmonella_ge.purity_filtered.txt", "BLUEPRINT_PE.T-cell_ge.purity_filtered.txt", "BLUEPRINT_SE.monocyte_ge.purity_filtered.txt", "BLUEPRINT_SE.neutrophil_ge.purity_filtered.txt"] #, "BrainSeq.brain_ge.purity_filtered.txt", "FUSION.adipose_naive_ge.purity_filtered.txt", "FUSION.muscle_naive_ge.purity_filtered.txt", "GENCORD.fibroblast_ge.purity_filtered.txt", "GENCORD.LCL_ge.purity_filtered.txt", "GENCORD.T-cell_ge.purity_filtered.txt", "GEUVADIS.LCL_ge.purity_filtered.txt"] #, "GTEx.adipose_subcutaneous_ge.purity_filtered.txt", "GTEx.adipose_visceral_ge.purity_filtered.txt", "GTEx.adrenal_gland_ge.purity_filtered.txt", "GTEx.artery_aorta_ge.purity_filtered.txt", "GTEx.artery_coronary_ge.purity_filtered.txt", "GTEx.artery_tibial_ge.purity_filtered.txt", "GTEx.blood_ge.purity_filtered.txt", "GTEx.brain_amygdala_ge.purity_filtered.txt", "GTEx.brain_anterior_cingulate_cortex_ge.purity_filtered.txt", "GTEx.brain_caudate_ge.purity_filtered.txt", "GTEx.brain_cerebellar_hemisphere_ge.purity_filtered.txt", "GTEx.brain_cerebellum_ge.purity_filtered.txt", "GTEx.brain_cortex_ge.purity_filtered.txt", "GTEx.brain_frontal_cortex_ge.purity_filtered.txt", "GTEx.brain_hippocampus_ge.purity_filtered.txt", "GTEx.brain_hypothalamus_ge.purity_filtered.txt", "GTEx.brain_nucleus_accumbens_ge.purity_filtered.txt", "GTEx.brain_putamen_ge.purity_filtered.txt", "GTEx.brain_spinal_cord_ge.purity_filtered.txt", "GTEx.brain_substantia_nigra_ge.purity_filtered.txt", "GTEx.breast_ge.purity_filtered.txt"] #, "GTEx.colon_sigmoid_ge.purity_filtered.txt", "GTEx.colon_transverse_ge.purity_filtered.txt", "GTEx.esophagus_gej_ge.purity_filtered.txt", "GTEx.esophagus_mucosa_ge.purity_filtered.txt", "GTEx.esophagus_muscularis_ge.purity_filtered.txt", "GTEx.fibroblast_ge.purity_filtered.txt", "GTEx.heart_atrial_appendage_ge.purity_filtered.txt", "GTEx.heart_left_ventricle_ge.purity_filtered.txt", "GTEx.LCL_ge.purity_filtered.txt", "GTEx.liver_ge.purity_filtered.txt", "GTEx.lung_ge.purity_filtered.txt", "GTEx.minor_salivary_gland_ge.purity_filtered.txt", "GTEx.muscle_ge.purity_filtered.txt", "GTEx.nerve_tibial_ge.purity_filtered.txt", "GTEx.ovary_ge.purity_filtered.txt", "GTEx.pancreas_ge.purity_filtered.txt", "GTEx.pituitary_ge.purity_filtered.txt", "GTEx.prostate_ge.purity_filtered.txt", "GTEx.skin_not_sun_exposed_ge.purity_filtered.txt", "GTEx.skin_sun_exposed_ge.purity_filtered.txt", "GTEx.small_intestine_ge.purity_filtered.txt", "GTEx.spleen_ge.purity_filtered.txt", "GTEx.stomach_ge.purity_filtered.txt", "GTEx.testis_ge.purity_filtered.txt", "GTEx.thyroid_ge.purity_filtered.txt", "GTEx.uterus_ge.purity_filtered.txt", "GTEx.vagina_ge.purity_filtered.txt", "HipSci.iPSC_ge.purity_filtered.txt", "Lepik_2017.blood_ge.purity_filtered.txt", "Nedelec_2016.macrophage_Listeria_ge.purity_filtered.txt", "Nedelec_2016.macrophage_naive_ge.purity_filtered.txt", "Nedelec_2016.macrophage_Salmonella_ge.purity_filtered.txt", "Quach_2016.monocyte_IAV_ge.purity_filtered.txt", "Quach_2016.monocyte_LPS_ge.purity_filtered.txt", "Quach_2016.monocyte_naive_ge.purity_filtered.txt", "Quach_2016.monocyte_Pam3CSK4_ge.purity_filtered.txt", "Quach_2016.monocyte_R848_ge.purity_filtered.txt", "ROSMAP.brain_naive_ge.purity_filtered.txt"] #, "Schmiedel_2018.B-cell_naive_ge.purity_filtered.txt", "Schmiedel_2018.CD4_T-cell_anti-CD3-CD28_ge.purity_filtered.txt", "Schmiedel_2018.CD4_T-cell_naive_ge.purity_filtered.txt", "Schmiedel_2018.CD8_T-cell_anti-CD3-CD28_ge.purity_filtered.txt", "Schmiedel_2018.CD8_T-cell_naive_ge.purity_filtered.txt", "Schmiedel_2018.monocyte_CD16_naive_ge.purity_filtered.txt", "Schmiedel_2018.monocyte_naive_ge.purity_filtered.txt", "Schmiedel_2018.NK-cell_naive_ge.purity_filtered.txt", "Schmiedel_2018.Tfh_memory_ge.purity_filtered.txt", "Schmiedel_2018.Th1-17_memory_ge.purity_filtered.txt", "Schmiedel_2018.Th17_memory_ge.purity_filtered.txt", "Schmiedel_2018.Th1_memory_ge.purity_filtered.txt", "Schmiedel_2018.Th2_memory_ge.purity_filtered.txt", "Schmiedel_2018.Treg_memory_ge.purity_filtered.txt", "Schmiedel_2018.Treg_naive_ge.purity_filtered.txt", "Schwartzentruber_2018.sensory_neuron_ge.purity_filtered.txt"] #, "TwinsUK.blood_ge.purity_filtered.txt", "TwinsUK.fat_ge.purity_filtered.txt", "TwinsUK.LCL_ge.purity_filtered.txt", "TwinsUK.skin_ge.purity_filtered.txt", "van_de_Bunt_2015.pancreatic_islet_ge.purity_filtered.txt"]
#ge_failinimede_list = ["HipSci.iPSC_ge.purity_filtered.txt", "Alasoo_2018.macrophage_IFNg_ge.purity_filtered.txt"]
#ge_failinimede_list = ["0ge.txt"]
ge_tabelite_pikkus = 0
ge_count_all = defaultdict(int)

# seda peab vaid 1 kord jooksutama
if not path.exists(r"C:\Users\pr\Documents\kool\6 semester 2020-21\Baka\andmed\\" + "ge_1.csv"):
    for fail in ge_failinimede_list: #algtöötleme kõik andmed läbi ning salvestame failidesse
        algtöötlus(fail)

#käime faili(kromosoomi)põhiselt kõik andmed läbi
for i in range(1): #TODO muidu 22
    ge_tabel = pd.read_csv(r"C:\Users\pr\Documents\kool\6 semester 2020-21\Baka\andmed\\" + "ge_" + str(i + 1) + ".csv")  # loeb tabeli sisse
    ge_tabelite_pikkus += len(ge_tabel)
    ge_count = ülekatete_leidmine(ge_tabel, 'ge', cs_suurus, False)

    for count in ge_count:
        ge_count_all[count] += ge_count[count]


print("Geeniekspressiooniga leitud ülekatted:", ge_count_all, ", algne tabeli suurus:", ge_tabelite_pikkus, ", leitud ülekatete protsent kogu andmetest:", (round((sum(ge_count_all.values()) / ge_tabelite_pikkus)*100, 5)), "%")



# RNA splaissimine (contained
print("\nRNA splaissimine:")
txrev_failinimede_list = ["Alasoo_2018.macrophage_IFNg+Salmonella_txrev.purity_filtered.txt", "Alasoo_2018.macrophage_IFNg_txrev.purity_filtered.txt", "Alasoo_2018.macrophage_naive_txrev.purity_filtered.txt", "Alasoo_2018.macrophage_Salmonella_txrev.purity_filtered.txt", "BLUEPRINT_PE.T-cell_txrev.purity_filtered.txt", "BLUEPRINT_SE.monocyte_txrev.purity_filtered.txt", "BLUEPRINT_SE.neutrophil_txrev.purity_filtered.txt"] #, "BrainSeq.brain_txrev.purity_filtered.txt", "FUSION.adipose_naive_txrev.purity_filtered.txt", "FUSION.muscle_naive_txrev.purity_filtered.txt", "GENCORD.fibroblast_txrev.purity_filtered.txt", "GENCORD.LCL_txrev.purity_filtered.txt", "GENCORD.T-cell_txrev.purity_filtered.txt", "GEUVADIS.LCL_txrev.purity_filtered.txt"] #, "GTEx.adipose_subcutaneous_txrev.purity_filtered.txt", "GTEx.adipose_visceral_txrev.purity_filtered.txt", "GTEx.adrenal_gland_txrev.purity_filtered.txt", "GTEx.artery_aorta_txrev.purity_filtered.txt", "GTEx.artery_coronary_txrev.purity_filtered.txt", "GTEx.artery_tibial_txrev.purity_filtered.txt", "GTEx.blood_txrev.purity_filtered.txt", "GTEx.brain_amygdala_txrev.purity_filtered.txt", "GTEx.brain_anterior_cingulate_cortex_txrev.purity_filtered.txt", "GTEx.brain_caudate_txrev.purity_filtered.txt", "GTEx.brain_cerebellar_hemisphere_txrev.purity_filtered.txt", "GTEx.brain_cerebellum_txrev.purity_filtered.txt", "GTEx.brain_cortex_txrev.purity_filtered.txt", "GTEx.brain_frontal_cortex_txrev.purity_filtered.txt", "GTEx.brain_hippocampus_txrev.purity_filtered.txt", "GTEx.brain_hypothalamus_txrev.purity_filtered.txt", "GTEx.brain_nucleus_accumbens_txrev.purity_filtered.txt", "GTEx.brain_putamen_txrev.purity_filtered.txt", "GTEx.brain_spinal_cord_txrev.purity_filtered.txt", "GTEx.brain_substantia_nigra_txrev.purity_filtered.txt", "GTEx.breast_txrev.purity_filtered.txt", "GTEx.colon_sigmoid_txrev.purity_filtered.txt", "GTEx.colon_transverse_txrev.purity_filtered.txt", "GTEx.esophagus_gej_txrev.purity_filtered.txt", "GTEx.esophagus_mucosa_txrev.purity_filtered.txt", "GTEx.esophagus_muscularis_txrev.purity_filtered.txt", "GTEx.fibroblast_txrev.purity_filtered.txt", "GTEx.heart_atrial_appendage_txrev.purity_filtered.txt", "GTEx.heart_left_ventricle_txrev.purity_filtered.txt", "GTEx.LCL_txrev.purity_filtered.txt", "GTEx.liver_txrev.purity_filtered.txt", "GTEx.lung_txrev.purity_filtered.txt", "GTEx.minor_salivary_gland_txrev.purity_filtered.txt", "GTEx.muscle_txrev.purity_filtered.txt", "GTEx.nerve_tibial_txrev.purity_filtered.txt", "GTEx.ovary_txrev.purity_filtered.txt", "GTEx.pancreas_txrev.purity_filtered.txt", "GTEx.pituitary_txrev.purity_filtered.txt", "GTEx.prostate_txrev.purity_filtered.txt", "GTEx.skin_not_sun_exposed_txrev.purity_filtered.txt", "GTEx.skin_sun_exposed_txrev.purity_filtered.txt", "GTEx.small_intestine_txrev.purity_filtered.txt", "GTEx.spleen_txrev.purity_filtered.txt", "GTEx.stomach_txrev.purity_filtered.txt", "GTEx.testis_txrev.purity_filtered.txt", "GTEx.thyroid_txrev.purity_filtered.txt", "GTEx.uterus_txrev.purity_filtered.txt", "GTEx.vagina_txrev.purity_filtered.txt", "HipSci.iPSC_txrev.purity_filtered.txt", "Lepik_2017.blood_txrev.purity_filtered.txt", "Nedelec_2016.macrophage_Listeria_txrev.purity_filtered.txt", "Nedelec_2016.macrophage_naive_txrev.purity_filtered.txt", "Nedelec_2016.macrophage_Salmonella_txrev.purity_filtered.txt", "Quach_2016.monocyte_IAV_txrev.purity_filtered.txt", "Quach_2016.monocyte_LPS_txrev.purity_filtered.txt", "Quach_2016.monocyte_naive_txrev.purity_filtered.txt", "Quach_2016.monocyte_Pam3CSK4_txrev.purity_filtered.txt", "Quach_2016.monocyte_R848_txrev.purity_filtered.txt", "ROSMAP.brain_naive_txrev.purity_filtered.txt"] #, "Schmiedel_2018.B-cell_naive_txrev.purity_filtered.txt", "Schmiedel_2018.CD4_T-cell_anti-CD3-CD28_txrev.purity_filtered.txt", "Schmiedel_2018.CD4_T-cell_naive_txrev.purity_filtered.txt", "Schmiedel_2018.CD8_T-cell_anti-CD3-CD28_txrev.purity_filtered.txt", "Schmiedel_2018.CD8_T-cell_naive_txrev.purity_filtered.txt", "Schmiedel_2018.monocyte_CD16_naive_txrev.purity_filtered.txt", "Schmiedel_2018.monocyte_naive_txrev.purity_filtered.txt", "Schmiedel_2018.NK-cell_naive_txrev.purity_filtered.txt", "Schmiedel_2018.Tfh_memory_txrev.purity_filtered.txt", "Schmiedel_2018.Th1-17_memory_txrev.purity_filtered.txt", "Schmiedel_2018.Th17_memory_txrev.purity_filtered.txt", "Schmiedel_2018.Th1_memory_txrev.purity_filtered.txt", "Schmiedel_2018.Th2_memory_txrev.purity_filtered.txt", "Schmiedel_2018.Treg_memory_txrev.purity_filtered.txt", "Schmiedel_2018.Treg_naive_txrev.purity_filtered.txt", "Schwartzentruber_2018.sensory_neuron_txrev.purity_filtered.txt"] #, "TwinsUK.blood_txrev.purity_filtered.txt", "TwinsUK.fat_txrev.purity_filtered.txt", "TwinsUK.LCL_txrev.purity_filtered.txt", "TwinsUK.skin_txrev.purity_filtered.txt", "van_de_Bunt_2015.pancreatic_islet_txrev.purity_filtered.txt"]
#txrev_failinimede_list = ["HipSci.iPSC_txrev.purity_filtered.txt", "Alasoo_2018.macrophage_IFNg_txrev.purity_filtered.txt"]
#txrev_failinimede_list = ["0txrev.txt"]
txrev_tabelite_pikkus = 0
txrev_tabelid = []
txrev_count_all = defaultdict(int)

# seda peab vaid 1 kord jooksutama

if not (path.exists(r"C:\Users\pr\Documents\kool\6 semester 2020-21\Baka\andmed\\" + "txrev_1.csv")):
    for fail in txrev_failinimede_list: #algtöötleme kõik andmed läbi ning salvestame failidesse
        algtöötlus(fail, "spalissimine")

#käime faili(kromosoomi)põhiselt kõik andmed läbi
for i in range(1): #TODO: 22

    txrev_tabel = pd.read_csv(r"C:\Users\pr\Documents\kool\6 semester 2020-21\Baka\andmed\\" + "txrev" + "_" + str(i + 1) + ".csv")  # loeb tabeli sisse
    txrev_tabelite_pikkus += len(txrev_tabel)
    txrev_count = ülekatete_leidmine(txrev_tabel, 'txrev', cs_suurus, True)

    for count in txrev_count:
        txrev_count_all[count] += txrev_count[count]


print("RNA splaissimisega leitud ülekatted:", txrev_count_all, ", algne tabeli suurus:", txrev_tabelite_pikkus, ", leitud ülekatete protsent kogu andmetest:", (round((sum(txrev_count_all.values()) / txrev_tabelite_pikkus)*100, 5)), "%")



#txrevise promootor (upstream)
#print("\nPromootor:")
#up_tabel = algtöötlus("HipSci.iPSC_txrev.purity_filtered.txt", "promootor")
#up_cs_list = cs_filtreerimine(up_tabel, cs_suurus, True)
#up_count = ülekatete_leidmine(up_tabel, up_cs_list, 'txrev_up')
#print("Promootoriga leitud ülekatted: ", up_count, ", algne tabeli suurus:", len(up_tabel), ", leitud ülekatete protsent kogu andmetest:", (round((sum(up_count.values()) / len(up_tabel))*100, 5)), "%")

ge_count_all, txrev_count_all, ge_labels, txrev_labels = ülekatete_korrastamine(ge_count_all, txrev_count_all, ge_tabelite_pikkus, txrev_tabelite_pikkus)


def joonista(ge, txrev, ge_labels, txrev_labels):
    #fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(16, 8))
    fig, (ax0) = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))

    # set width of bar
    barWidth = 0.45

    # set height of bar
    bars1 = list(ge.values())
    bars2 = list(txrev.values())

    # Set position of bar on X axis
    r1 = np.arange(len(bars1))
    r2 = [x + barWidth for x in r1]

    # Make the plot
    ax0.bar(r1, bars1, color='c', width=barWidth, edgecolor='white', label='Geeniekspressioon')
    ax0.bar(r2, bars2, color='coral', width=barWidth, edgecolor='white', label='RNA splaissimine')

    ax0.set_title('Tulemuste võrdlus, cs <' + str(cs_suurus))

    # Add xticks on the middle of the group bars
    ax0.set_xlabel('Ülekatted')
    ax0.set_ylabel('Geenide arv protsentuaalselt kogu andmestikust')

    ax0.set_xticks([r + (barWidth / 2) for r in range(len(bars1))])
    ax0.set_xticklabels(list(ge.keys()))
    i = 0
    labels = ge_labels + txrev_labels

    for p in ax0.patches:
        ax0.annotate(str(labels[i]), (p.get_x() * 1.005 + barWidth / 6, p.get_height() * 1.005))
        i += 1

    ax0.legend()

    # ge.pop(1, None)
    # txrev.pop(1, None)
    #
    # # set height of bar
    # bars1 = list(ge.values())
    # bars2 = list(txrev.values())
    #
    # # Set position of bar on X axis
    # r1 = np.arange(len(bars1))
    # r2 = [x + barWidth for x in r1]
    #
    # # Make the plot
    # ax1.bar(r1, bars1, color='c', width=barWidth, edgecolor='white', label='Geeniekspressioon')
    # ax1.bar(r2, bars2, color='coral', width=barWidth, edgecolor='white', label='RNA splaissimine')
    #
    # ax1.set_title('Tulemuste võrdlus 2+ ülekattega')
    #
    # # Add xticks on the middle of the group bars
    # ax1.set_xlabel('Ülekatteid')
    # ax1.set_ylabel('Geenide arv')
    # ax1.set_xticks([r + (barWidth / 2) for r in range(len(bars1))])
    # ax1.set_xticklabels(list(ge.keys()))
    #
    # # Annotate every column
    # i = 0
    # labels = ge_labels[1:] + txrev_labels[1:] #ühe ülekattega veerud jätame välja
    # for p in ax1.patches:
    #     ax1.annotate(str(labels[i]), (p.get_x() * 1.005 + barWidth / 6, p.get_height() * 1.005))
    #     i += 1
    # ax1.legend()


    # CShow graphic
    plt.show()


joonista(ge_count_all, txrev_count_all, ge_labels, txrev_labels)
