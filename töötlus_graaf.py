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
def algtöötlus(failinimi, salvesta, vaatleme_kromosoome, tüüp):
    tabel = pd.read_csv(r"C:\Users\pr\Documents\kool\6 semester 2020-21\Baka\andmed\eraldi\\" + failinimi, sep="	")  # loeb faili sisse
    tabel = tabel[['phenotype_id', 'chr', 'cs_id', 'pip', 'cs_size', 'variant_id']]  # mis veerud jätta alles
    tabel = tabel.rename(columns={'chr': 'Chromosome'})
    tabel = pd.concat([tabel, tabel["variant_id"].str.split("_", expand=True)[1].rename("Start")], axis=1)  # uus veerg Start
    tabel["Start"] = pd.to_numeric(tabel["Start"])
    tabel['End'] = tabel['Start'] + 1  # Pyranges puhul peab mingi algus ja lõppkoht olema
    tabel = tabel[['phenotype_id', 'variant_id', 'cs_id', 'Chromosome', 'Start', 'End', 'pip', 'cs_size']]  # mis veerud jätta alles

    if tüüp == 2:  # vaid txrev puhul
        tabel = tabel[tabel['phenotype_id'].str.contains("contained")]  # jätab alles vaid andmed splaissimise kohta

    elif tüüp == 3:  # vaid txrev puhul promootor
        tabel = tabel[tabel['phenotype_id'].str.contains("upstream")]

    elif tüüp == 4:  # vaid txrev puhul downstream
        tabel = tabel[tabel['phenotype_id'].str.contains("downstream")]

    tabel["cs_id"] = tabel["cs_id"] + "_" + failinimi.split(".purity_filtered.txt")[0] #lisame cs_id veerule andmebaasi nime juurde

    #teeme igale kromosoomile oma tabeli
    tabelid = list(range(vaatleme_kromosoome)) #TODO muidu 22

    for i in range(vaatleme_kromosoome): #kirjutame iga kromosoomi oma faili #TODO muidu 22
        tabelid[i] = (tabel[(tabel['Chromosome'] == 1 + i)])
        if path.exists(r"C:\Users\pr\Documents\kool\6 semester 2020-21\Baka\andmed\koos\\" + salvesta + str(i + 1) + ".csv"):
            tabelid[i].to_csv(
                r"C:\Users\pr\Documents\kool\6 semester 2020-21\Baka\andmed\koos\\" + salvesta + str(i + 1) + ".csv", mode='a', index=False, header=False)
        else:
            tabelid[i].to_csv(
                r"C:\Users\pr\Documents\kool\6 semester 2020-21\Baka\andmed\koos\\" + salvesta + str(i + 1) + ".csv", mode='a', index=False)  # esimene kord salvestame koos veerunimetustega


'''
Tabelist ülekatete leidmine kasutades Pyranges klasterdamist. Joonistab klastrid, kus tippudeks on cs-id ja servad on nende vahel siis kui need jagavad omavahel samu variante
Sisend: tabel, faili nime algus, maksimaalne credible set-i suurus, kas on tegemist splaissimisega
Väljund: leitud ülekatted
'''
def ülekatete_leidmine(tabel, tüüp, laiend, cs_suurus, kromosoom):
    tabel = tabel[tabel['cs_size'] <= cs_suurus]  # credible set maksimaalne suurus
    if tüüp == 2 or tüüp == 3 or tüüp == 4:
        tabel['phenotype_id'] = tabel['phenotype_id'].str.split(".", expand=True)  # jätab alles vaid geeni nime

    pr_tabel = pr.PyRanges(tabel)  # teeb PyRanges tabeli
    cluster = pr_tabel.cluster(by="Start") #klasterdame asukohapõhiselt
    cluster_df = cluster.df

    # joined = pr_tabel.join(pr_tabel, report_overlap=True)
    # joined_df = joined.df
    #print(joined_df)

    if (kromosoom == 1):
        # print(cluster_df)
        cluster_df.to_csv(r"tulemused\\" + laiend + "cluster_tulemus.csv")

        merge = pr_tabel.merge(count=True, by="Start")
        merge_df = merge.df
        merge_df.to_csv(r"tulemused\\" + laiend + "merge_tulemus.csv")

    # cluster_df = cluster_df[:50] #TODO: eemalda see

    #saame kätte cs-id
    cs_id = set((tabel.cs_id.values))

    g = ig.Graph()

    # tippude lisamine
    g.add_vertices(len(cs_id))
    eelmine_klaster = -1
    # lisa_servad = []
    i = 0
    cs_id_indeks = {}
    for cs in cs_id:
        cs_id_indeks[cs] = i
        i += 1

    eelmise_indeks = -1

    #print(cs_id_indeks)
    for i, rida in cluster_df.iterrows():
        klaster = rida.Cluster
        #lisame tipu
        #print("tipp", i)
        #print(rida.phenotype_id, " ", rida.cs_id, " ", rida.Cluster)

        indeks = cs_id_indeks[rida.cs_id]

        try:
            g.vs[indeks]["positions"].append(rida.Start)

        except:
            g.vs[indeks]['cs'] = rida.cs_id
            g.vs[indeks]["positions"] = [rida.Start]
            g.vs[indeks]["phenotype_id"] = rida.phenotype_id


        #lisame servad
        # lisame servad vaid eelmise vahel, kuna lõpuks see tulemust ei muuda (muuda arvutuskäiku aeglasemaks)
        if (klaster == eelmine_klaster):
            g.add_edge(indeks, eelmise_indeks)
        else:
            eelmine_klaster = klaster
        eelmise_indeks = indeks
        # if (klaster == eelmine_klaster): #kuulub eelnevaga samasse klastrisse
        #     lisa_servad.append(indeks)
        # else:
        #     #teeme servad ühe klastri vahel ära
        #     for servad in combinations(lisa_servad, 2):
        #         g.add_edge(servad[0], servad[1])
        #
        #     #muudame praeguse klastri informatsiooni
        #     lisa_servad = [indeks]
        #     eelmine_klaster = klaster


    #kui lõpus jäi mõni serv lisamata
    if (klaster == eelmine_klaster):
        g.add_edge(indeks, eelmise_indeks)

    # for servad in combinations(lisa_servad, 2):
    #     g.add_edge(servad[0], servad[1])

    # moodustame klikke. Ühte klikki kuuluvad sellised tipud, mille iga tipu vahel leidub tee
    components = g.clusters()
    geenid_ülekattes = []
    komponentide_suurused = defaultdict(int)
    komponentide_geenid = []

    #leiame, mis geenid kuuluvad samadesse põhjuslike variantide klikkidesse (credible component)
    for klikk in components:
        geenid = set()
        #leiame geenid komponentidest
        for tipp in klikk:
            geenid.add(g.vs[tipp]["phenotype_id"])

        komponentide_suurused[len(klikk)] += 1
        komponentide_geenid.append(geenid)
        geenid_ülekattes.append(geenid)

    print(kromosoom)


    #salvestame leitud ülekatted
    ülekatteid = defaultdict(int)
    for komponent in geenid_ülekattes:
        ülekatteid[len(komponent)] += 1


    # if (kromosoom == 1):
    #     #salvestame tipud faili
    #     tipud_df = g.get_vertex_dataframe()
    #     #print(tipud_df)
    #     tipud_df.to_csv(r"tulemused\\" + nimi + "_tipud.csv")
    #
    #     #joonistame
    #     #g_subgraph = g
    #     g_subgraph = g.subgraph(list(range(100)))
    #     components = g_subgraph.clusters()
    #     layout = g_subgraph.layout(layout='auto')
    #     ig.plot(components)

    return ülekatteid, komponentide_suurused, komponentide_geenid


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



def main(tüüp, failinimede_list):
    vaatleme_kromosoome = 22
    tabelite_pikkus = 0
    ülekatteid_kokku = defaultdict(int)
    komponentide_suurused_kokku = {}
    komponentidesse_kuuluvad_geenid_kokku = []

    if tüüp == 1: #geeniekspressioon
        laiend = "ge"
        salvesta = "ge_"
        nimetus = "Geeniekspressiooniga"
    elif tüüp == 2: #RNA splaissimine
        laiend = "txrev"
        salvesta = "txrev_"
        nimetus = "RNA splaissimisega"
    elif tüüp == 3: #promootor
        laiend = "txrev"
        salvesta = "upstream_"
        nimetus = "Upstream-iga"
    elif tüüp == 4:
        laiend = "txrev"
        salvesta = "downstream_"
        nimetus = "Downstream-iga"
    else:
        laiend = "muu"
        salvesta = "muu_"
        nimetus = "Muuga"

    # seda peab vaid 1 kord jooksutama
    if not path.exists(r"C:\Users\pr\Documents\kool\6 semester 2020-21\Baka\andmed\koos\\" + salvesta + "1.csv"): #pole läbi jooksutatud
        for fail in failinimede_list:  # algtöötleme kõik andmed läbi ning salvestame failidesse
            algtöötlus(fail + laiend + ".purity_filtered.txt", salvesta, vaatleme_kromosoome, tüüp)


    for i in range(vaatleme_kromosoome):
        tabel = pd.read_csv(r"C:\Users\pr\Documents\kool\6 semester 2020-21\Baka\andmed\koos\\" + salvesta + str(i + 1) + ".csv")  # loeb tabeli sisse
        tabelite_pikkus += len(tabel)
        ülekatteid, komponentide_suurused, komponentidesse_kuuluvad_geenid = ülekatete_leidmine(tabel, tüüp, salvesta, cs_suurus, i + 1)

        komponentide_suurused_kokku.update(dict(komponentide_suurused))
        komponentidesse_kuuluvad_geenid_kokku.append(komponentidesse_kuuluvad_geenid)

        for ülekate in ülekatteid.keys():
            ülekatteid_kokku[ülekate] += ülekatteid[ülekate]


    mitu_geeni = 0
    for ülekate in ülekatteid_kokku.keys():
        if (ülekate > 1):
            mitu_geeni += ülekatteid_kokku[ülekate]


    # komponentide_suurused
    with open(file=r"tulemused\\" + laiend + "komponentide_suurused.txt", mode='w') as f:
        f.write(str(komponentide_suurused_kokku))
    with open(file=r"tulemused\\" + laiend + "komponentidesse_kuuluvad_geenid.txt", mode='w') as f:
        f.write(str(komponentidesse_kuuluvad_geenid_kokku))



    print("\n"+ nimetus + " leitud ülekatted:", dict(ülekatteid_kokku), ", algne tabeli suurus:",
          tabelite_pikkus, )  # ", leitud ülekatete protsent kogu andmetest:", (round((sum(ge_count_all.values()) / ge_tabelite_pikkus)*100, 5)), "%")
    print("Ühes klikis geene: üks -", ülekatteid_kokku[1], ", mitu -", mitu_geeni)
    print("Leitud geenide suhe: ", (round((mitu_geeni / ülekatteid_kokku[1]) * 100, 3)), "%")
    #print("Komponentide suurused:", komponentide_suurused_kokku)
    #print("Ühte komponenti kuuluvad geenid:", komponentidesse_kuuluvad_geenid_kokku)
    ülekatteid_üldine = {1: ülekatteid_kokku[1], 2: mitu_geeni}

    return ülekatteid_üldine


cs_suurus = 10  # credible set maksimaalne suurus
failinimede_list = ["Alasoo_2018.macrophage_IFNg+Salmonella_", "Alasoo_2018.macrophage_IFNg_", "Alasoo_2018.macrophage_naive_", 'BLUEPRINT_PE.T-cell_', 'BLUEPRINT_SE.monocyte_', 'BLUEPRINT_SE.neutrophil_', 'BrainSeq.brain_', 'FUSION.adipose_naive_', 'FUSION.muscle_naive_', 'GENCORD.fibroblast_', 'GENCORD.LCL_', 'GENCORD.T-cell_', 'GEUVADIS.LCL_', 'GTEx.adipose_subcutaneous_', 'GTEx.adipose_visceral_', 'GTEx.adrenal_gland_', 'GTEx.artery_aorta_', 'GTEx.artery_coronary_', 'GTEx.artery_tibial_', 'GTEx.blood_', 'GTEx.brain_amygdala_', 'GTEx.brain_anterior_cingulate_cortex_', 'GTEx.brain_caudate_', 'GTEx.brain_cerebellar_hemisphere_', 'GTEx.brain_cerebellum_', 'GTEx.brain_cortex_', 'GTEx.brain_frontal_cortex_', 'GTEx.brain_hippocampus_', 'GTEx.brain_hypothalamus_', 'GTEx.brain_nucleus_accumbens_', 'GTEx.brain_putamen_', 'GTEx.brain_spinal_cord_', 'GTEx.brain_substantia_nigra_', 'GTEx.breast_', 'GTEx.colon_sigmoid_', 'GTEx.colon_transverse_', 'GTEx.esophagus_gej_', 'GTEx.esophagus_mucosa_', 'GTEx.esophagus_muscularis_', 'GTEx.fibroblast_', 'GTEx.heart_atrial_appendage_', 'GTEx.heart_left_ventricle_', 'GTEx.LCL_', 'GTEx.liver_', 'GTEx.lung_', 'GTEx.minor_salivary_gland_', 'GTEx.muscle_', 'GTEx.nerve_tibial_', 'GTEx.ovary_', 'GTEx.pancreas_', 'GTEx.pituitary_', 'GTEx.prostate_', 'GTEx.skin_not_sun_exposed_', 'GTEx.skin_sun_exposed_', 'GTEx.small_intestine_', 'GTEx.spleen_', 'GTEx.stomach_', 'GTEx.testis_', 'GTEx.thyroid_', 'GTEx.uterus_', 'GTEx.vagina_', 'HipSci.iPSC_', 'Lepik_2017.blood_', 'Nedelec_2016.macrophage_Listeria_', 'Nedelec_2016.macrophage_naive_', 'Nedelec_2016.macrophage_Salmonella_', 'Quach_2016.monocyte_IAV_', 'Quach_2016.monocyte_LPS_', 'Quach_2016.monocyte_naive_', 'Quach_2016.monocyte_Pam3CSK4_', 'Quach_2016.monocyte_R848_', 'ROSMAP.brain_naive_', 'Schmiedel_2018.B-cell_naive_', 'Schmiedel_2018.CD4_T-cell_anti-CD3-CD28_', 'Schmiedel_2018.CD4_T-cell_naive_', 'Schmiedel_2018.CD8_T-cell_anti-CD3-CD28_', 'Schmiedel_2018.CD8_T-cell_naive_', 'Schmiedel_2018.monocyte_CD16_naive_', 'Schmiedel_2018.monocyte_naive_', 'Schmiedel_2018.NK-cell_naive_', 'Schmiedel_2018.Tfh_memory_', 'Schmiedel_2018.Th1-17_memory_', 'Schmiedel_2018.Th17_memory_', 'Schmiedel_2018.Th1_memory_', 'Schmiedel_2018.Th2_memory_', 'Schmiedel_2018.Treg_memory_', 'Schmiedel_2018.Treg_naive_', 'Schwartzentruber_2018.sensory_neuron_', 'TwinsUK.blood_', 'TwinsUK.fat_', 'TwinsUK.LCL_', 'TwinsUK.skin_', 'van_de_Bunt_2015.pancreatic_islet_']
#failinimede_list = ["HipSci.iPSC_", "Alasoo_2018.macrophage_IFNg+Salmonella_", "Alasoo_2018.macrophage_IFNg_", "Alasoo_2018.macrophage_naive_"]

# geeniekspressioon
print("\nGeeniekspressioon:")
ge_ülekatteid = main(1, failinimede_list)


# RNA splaissimine (contained
print("\nRNA splaissimine:")
txrev_ülekatteid = main(2, failinimede_list)


# txrevise promootor (upstream)
print("\nUpstream:")
upstream_ülekatteid = main(3, failinimede_list)

# txrevise  (downstream)
print("\nDownstream:")
downstream_ülekatteid = main(4, failinimede_list)


def joonista_kõik(ge, txrev, ge_labels, txrev_labels):
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

# ge_count_all, txrev_count_all, ge_labels, txrev_labels = ülekatete_korrastamine(ge_count_all, txrev_count_all, ge_tabelite_pikkus, txrev_tabelite_pikkus)
# joonista_kõik(ge_count_all, txrev_count_all, ge_labels, txrev_labels)


def joonista_üldine(ge, txrev, upstream, downstream):
    fig, (ax0) = plt.subplots( figsize=(8, 8))

    # set width of bar
    barWidth = 0.2

    # set height of bar
    bars1 = list(ge.values())
    bars2 = list(txrev.values())
    bars3 = list(upstream.values())
    bars4 = list(downstream.values())


    # Set position of bar on X axis
    r1 = np.arange(len(bars1))
    r2 = [x + barWidth for x in r1]
    r3 = [x + barWidth for x in r2]
    r4 = [x + barWidth for x in r3]

    # Make the plot
    ax0.bar(r1, bars1, color='c', width=barWidth, edgecolor='white', label='Geeniekspressioon')
    ax0.bar(r2, bars2, color='salmon', width=barWidth, edgecolor='white', label='RNA splaissimine')
    ax0.bar(r3, bars3, color='orchid', width=barWidth, edgecolor='white', label='Upstream')
    ax0.bar(r4, bars4, color='lightgreen', width=barWidth, edgecolor='white', label='Downstream')

    ax0.set_title('Geenide arv klikis vastavalt uurimistüübile, cs <=' + str(cs_suurus))

    # Add xticks on the middle of the group bars
    ax0.set_xlabel('Geenide arv klikis')
    ax0.set_ylabel('Klikkide arv')

    ax0.set_xticks([r + (barWidth / 4) for r in range(len(bars1))])
    ax0.set_xticklabels(["üks", "mitu"])
    i = 0
    labels = list(ge.values()) + list(txrev.values()) + list(upstream.values()) + list(downstream.values())

    for p in ax0.patches:
        ax0.annotate(str(labels[i]), (p.get_x() * 1.005 + barWidth / 6, p.get_height() * 1.005))
        i += 1

    ax0.legend()
    plt.show()


joonista_üldine(ge_ülekatteid, txrev_ülekatteid, upstream_ülekatteid, downstream_ülekatteid)
