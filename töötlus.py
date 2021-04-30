import pandas as pd
import pyranges as pr
import numpy as np
import matplotlib.pyplot as plt
import collections
from collections import defaultdict
from os import path
import igraph as ig
from itertools import combinations

pd.set_option('display.max_columns', None)
desired_width = 1000
pd.set_option('display.width', desired_width)
pd.set_option('display.max_colwidth', 100)
np.set_printoptions(linewidth=desired_width)
pd.options.mode.chained_assignment = None


def algtöötlus(failinimi, salvesta, vaatleme_kromosoome, tüüp):
    """
    Failist andmete sisselugemine, veergude lisamine ja kustutamine ning kromosoomipõhiselt failidesse salvestamine.
    :param failinimi: vaadeldava faili nimetus
    :param salvesta: faililaiend salvestamiseks
    :param vaatleme_kromosoome: mitmenda kromosoomini vaadeldakse
    :param tüüp: vaadeldav protsess (geeniekspressioon, transkriptsioon, splaissimine, promootor, terminaator)
    """
    tabel = pd.read_csv(r"C:\Users\pr\Documents\kool\6 semester 2020-21\Baka\andmed\eraldi\\" + failinimi, sep="	")  # loeb faili sisse
    tabel = tabel[['phenotype_id', 'chr', 'cs_id', 'pip', 'cs_size', 'variant_id']]  # mis veerud jätta alles
    tabel = tabel.rename(columns={'chr': 'Chromosome'})
    tabel = pd.concat([tabel, tabel["variant_id"].str.split("_", expand=True)[1].rename("Start")], axis=1)  # uus veerg Start
    tabel["Start"] = pd.to_numeric(tabel["Start"])
    tabel['End'] = tabel['Start'] + 1  # Pyranges puhul peab mingi algus ja lõppkoht olema
    tabel = tabel[['phenotype_id', 'variant_id', 'cs_id', 'Chromosome', 'Start', 'End', 'pip', 'cs_size']]  # mis veerud jätta alles

    # txrev faili puhul tuleb splaissimine, promootori ja terminaatori kasutuse andmed üksteisest eraldada
    if tüüp == 2:
        tabel = tabel[tabel['phenotype_id'].str.contains("contained")]

    elif tüüp == 3:
        tabel = tabel[tabel['phenotype_id'].str.contains("upstream")]

    elif tüüp == 4:
        tabel = tabel[tabel['phenotype_id'].str.contains("downstream")]

    tabel["cs_id"] = tabel["cs_id"] + "_" + failinimi.split(".purity_filtered.txt")[0]  # lisame cs_id veerule andmebaasi nime juurde

    # teeme igale kromosoomile oma tabeli
    tabelid = list(range(vaatleme_kromosoome))

    for i in range(vaatleme_kromosoome):  # kirjutame iga kromosoomi oma faili
        tabelid[i] = (tabel[(tabel['Chromosome'] == 1 + i)])
        if path.exists(
                r"C:\Users\pr\Documents\kool\6 semester 2020-21\Baka\andmed\koos\\" + salvesta + str(i + 1) + ".csv"):
            tabelid[i].to_csv(
                r"C:\Users\pr\Documents\kool\6 semester 2020-21\Baka\andmed\koos\\" + salvesta + str(i + 1) + ".csv",
                mode='a', index=False, header=False)
        else:
            tabelid[i].to_csv(
                r"C:\Users\pr\Documents\kool\6 semester 2020-21\Baka\andmed\koos\\" + salvesta + str(i + 1) + ".csv",
                mode='a', index=False)  # esimene kord salvestame koos veerunimetustega


def ülekatete_leidmine(tabel, tüüp, laiend, cs_suurus, kromosoom, transkript_sõnastik):
    """
    Tabelist ülekatete leidmine kasutades Pyranges klasterdamist. Koostab ühendatud komponente, kus tippudeks on CSid ja servad on nende vahel siis kui need jagavad omavahel samu variante
    :param tabel: mis andeid analüüsitakse
    :param tüüp: vaadeldav protsess (geeniekspressioon, transkriptsioon, splaissimine, promootor, terminaator)
    :param laiend: mis laiendiga faili sisse lugeda
    :param cs_suurus: maksimaalne CS suurus
    :param kromosoom: mitmendat kromosoomi vaatab
    :param transkript_sõnastik: kui tüüp = transkript, siis sisaldab sõnastikku, et saada kätte geeni_id, muidu tühi
    :return: sõnastik, mitu geeni iga komponent sisaldab
    """
    tabel = tabel[tabel['cs_size'] <= cs_suurus]  # credible set maksimaalne suurus

    if tüüp == 2 or tüüp == 3 or tüüp == 4:
        tabel['phenotype_id'] = tabel['phenotype_id'].str.split(".", expand=True)  # jätab alles vaid geeni_id

    if tüüp == 5:
        tabel["gene_id"] = tabel.apply(lambda row: transkript_sõnastik[row.phenotype_id], axis=1)  # saame igale phenotype_id-le vastava geeni_id
    else:
        tabel["gene_id"] = tabel.phenotype_id  # muudel juhtudel on geeni_id teada
    tabel.drop(["phenotype_id"], axis=1)

    pr_tabel = pr.PyRanges(tabel)  # teeb PyRanges tabeli
    cluster = pr_tabel.cluster(by="Start")  # klasterdab asukohapõhiselt
    cluster_df = cluster.df  # teeb Pandas DataFrame'i


    cluster_df.to_csv(r"tulemused\\" + str(cs_suurus) + "\\" + laiend + "cluster_tulemus.csv", mode='a', index=False)  # kirjutab faili

    # annab iga variandi kohta teada, mitmes CSis esineb
    # merge = pr_tabel.merge(count=True, by="Start")
    # merge_df = merge.df
    # merge_df.to_csv(r"tulemused\\" + str(cs_suurus) + "\\" + laiend + "merge_tulemus.csv", mode='a', index=False) # kirjutab faili

    # saab kätte unikaalsed CSide nimetused
    cs_id = set((tabel.cs_id.values))

    g = ig.Graph()  # loob graafi

    # tippude loomine
    g.add_vertices(len(cs_id))
    eelmine_klaster = -1
    i = 0
    cs_id_indeks = {}  # tippude (CSide) indeksite salvestamine
    for cs in cs_id:
        cs_id_indeks[cs] = i
        i += 1

    eelmise_indeks = -1

    for i, rida in cluster_df.iterrows():  # vaatab klasterdamisest saadud tabeli reahaaval läbi
        klaster = rida.Cluster
        indeks = cs_id_indeks[rida.cs_id]

        try:  # kui selle tipu (CSi) kohta on juba tipp loodud
            g.vs[indeks]["positions"].append(rida.Start)  # lisab veel ühe variandi asukoha, mis kuulub sellesse CSi

        except:  # muidu lisatakse tipp ja selle info
            g.vs[indeks]['cs'] = rida.cs_id
            g.vs[indeks]["positions"] = [rida.Start]
            g.vs[indeks]["gene_id"] = rida.gene_id

        # servade lisamine
        # servad lisatakse ühes klastris vaid kahe järjestikuse CSi vahel (ei muuda lõpus tulemust ja muudab arvutuskäigu kiiremaks
        if klaster == eelmine_klaster:
            g.add_edge(indeks, eelmise_indeks)
        else:
            eelmine_klaster = klaster
        eelmise_indeks = indeks

        # kui ikkagi vaja saada servad kõikide CSide vahel kätte
    #     lisa_servad = [] # milliste CSide vahel teha servi
    #     if (klaster == eelmine_klaster): # kuulub eelnevaga samasse klastrisse
    #         lisa_servad.append(indeks)
    #     else:
    #         # vaadeldav on uus klaster, teeb eelmise klastri servad ära
    #         for servad in combinations(lisa_servad, 2):
    #             g.add_edge(servad[0], servad[1])
    #
    #         # muudab praeguse klastri informatsiooni
    #         lisa_servad = [indeks] # vaadeldava CSi indeks salvestatakse listi
    #         eelmine_klaster = klaster
    #
    # # joonistab servad viimase klastri kohta ära
    # for servad in combinations(lisa_servad, 2):
    #     g.add_edge(servad[0], servad[1])

    komponendid = g.clusters()  # leiab ühendatud komponendid (sidusad alamgraafid)
    komponentide_geenid = []  # millised geenid kuuluvad samasse ühendatud komponenti

    # leiab, mis geenid kuuluvad samadesse ühendatud komponentidesse
    for komponent in komponendid:
        geenid = set()
        # leiab kõik komponendi geenid
        for tipp in komponent:
            geenid.add(g.vs[tipp]["gene_id"])

        komponentide_geenid.append(geenid)  # mis geenid kuuluvad samasse komponenti

    # salvestab leitud ülekatted
    ülekatted = defaultdict(int)
    for komponent in komponentide_geenid:
        ülekatted[len(komponent)] += 1

    # andmete faili salvestamine
    with open(file=r"tulemused\\" + str(cs_suurus) + "\\" + laiend + "komponentidesse_kuuluvad_geenid.txt",
              mode='a') as f:
        f.write(str(kromosoom) + ": " + str(komponentide_geenid) + "\n")

    return ülekatted


def main(tüüp, failinimede_list):
    """
    Leiab, kui palju on vaadeldava protsessi puhul üht või mitut geeni sisaldavaid ühendatud komponente.
    Prindib huvipakkuva informatsiooni välja.
    :param tüüp: vaadeldav protsess (geeniekspressioon, transkriptsioon, splaissimine, promootor, terminaator)
    :param failinimede_list: mis failid vaadeldakse läbi
    :return: kui palju oli vaid üht geeni sisaldavad ühendatud komponente ja kui palju oli rohkem kui üht geeni sisaldavaid komponente
    """
    vaatleme_kromosoome = 22
    ülekatteid_kokku = defaultdict(int)  # kui suuri ülekatteid on leitud

    if tüüp == 1:  # geeniekspressioon
        laiend = "ge"
        salvesta = "ge_"
        nimetus = "Geeniekspressiooniga"
    elif tüüp == 2:  # RNA splaissimine
        laiend = "txrev"
        salvesta = "txrev_"
        nimetus = "RNA splaissimisega"
    elif tüüp == 3:  # promootor
        laiend = "txrev"
        salvesta = "upstream_"
        nimetus = "Promootori kasutusega "
    elif tüüp == 4:  # terminaator
        laiend = "txrev"
        salvesta = "downstream_"
        nimetus = "Terminaatori kasutusega"
    elif tüüp == 5:  # transkriptsioon
        laiend = "tx"
        salvesta = "tx_"
        nimetus = "Transkriptidega"
    else:
        laiend, salvesta, nimetus = "muu"

    # seda peab vaid 1 kord jooksutama, et salvestada kogu informatsioon kromosoomipõhiselt
    if not path.exists(
            r"C:\Users\pr\Documents\kool\6 semester 2020-21\Baka\andmed\koos\\" + salvesta + "1.csv"):  # pole läbi jooksutatud
        for fail in failinimede_list:  # algtöötleme kõik andmed läbi ning salvestame failidesse
            algtöötlus(fail + laiend + ".purity_filtered.txt", salvesta, vaatleme_kromosoome, tüüp)

    if tüüp == 5:
        # vajalik transcriptide fenotüüpide teisendamiseks geenideks
        metadata = pd.read_csv(
            r"C:\Users\pr\Documents\kool\6 semester 2020-21\Baka\andmed\eraldi\transcript_usage_Ensembl_96_phenotype_metadata.tsv",
            sep='\t', header=0)  # loeb metadata sisse
        metadata = metadata[["phenotype_id", "gene_id"]]
        transkript_sõnastik = metadata.set_index("phenotype_id").to_dict()["gene_id"]
    else:
        transkript_sõnastik = {}

    # analüüsib kromosoomipõhiselt kõik läbi
    for i in range(vaatleme_kromosoome):
        tabel = pd.read_csv(r"C:\Users\pr\Documents\kool\6 semester 2020-21\Baka\andmed\koos\\" + salvesta + str(i + 1) + ".csv")  # loeb tabeli sisse
        ülekatteid = ülekatete_leidmine(tabel, tüüp, salvesta, cs_suurus, i + 1, transkript_sõnastik)  # leiab ülekatted

        for ülekate in ülekatteid.keys():
            ülekatteid_kokku[ülekate] += ülekatteid[ülekate]  # salvestab, palju on ülekatteid kokku leitud

    mitu_geeni = 0  # palju on rohkem kui üht geeni sisaldavaid komponente
    for ülekate in ülekatteid_kokku.keys():
        if ülekate > 1:
            mitu_geeni += ülekatteid_kokku[ülekate]

    ülekatteid_kokku = dict(collections.OrderedDict(sorted(ülekatteid_kokku.items())))  # järjestab sõnastiku

    komponente_kokku = 0  # komponentide arv kokku
    paarid = [(k, ülekatteid_kokku[k]) for k in ülekatteid_kokku]  # paarid: geenide arv komponendis ja kui palju selliseid komponente leidus
    for paar in paarid:  # leiab komponentide arvu
        komponente_kokku += paar[0] * paar[1]

    print(nimetus + " leitud ülekatted:", dict(ülekatteid_kokku))
    print("Ühes ühendatud komponendis geene: üks -", ülekatteid_kokku[1], " mitu -", mitu_geeni)
    print("Vaid üht geeni sisaldavate ühendatud komponentide suhe kogu komponentide arvuga: ",
          (round((ülekatteid_kokku[1] / (ülekatteid_kokku[1] + mitu_geeni)) * 100, 3)), "%")
    print("Ühendatud komponente kokku ", komponente_kokku)
    ülekatteid_üldine = {1: ülekatteid_kokku[1], 2: mitu_geeni}

    return ülekatteid_üldine


# kõik töös kasutatud andmed
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
# failinimede_list = ["HipSci.iPSC_", "Alasoo_2018.macrophage_IFNg+Salmonella_"]  # katsetamiseks


ülekatted = []
for cs in range(10, 101, 10):  # mis CSi vahemikus koodi läbi jooksutada
    cs_suurus = cs
    print("\n\nCS suurus:", cs_suurus)

    # ge
    print("\nGeeniekspressioon:")
    ge_ülekatteid = main(1, failinimede_list)

    # txrevise contained
    print("\nRNA splaissimine:")
    txrev_ülekatteid = main(2, failinimede_list)

    # txrevise upstream
    print("\nPromootor:")
    upstream_ülekatteid = main(3, failinimede_list)

    # txrevise downstream
    print("\nTerminaator:")
    downstream_ülekatteid = main(4, failinimede_list)

    # tx
    print("\nTranskript:")
    tx_ülekatteid = main(5, failinimede_list)


def joonista(ge, txrev, upstream, downstream, transcript):
    fig, (ax0) = plt.subplots(figsize=(8, 8))

    # tulba laius
    bar_width = 0.18

    # tulba kõrgus
    bars1 = list(ge.values())
    bars2 = list(txrev.values())
    bars3 = list(upstream.values())
    bars4 = list(downstream.values())
    bars5 = list(transcript.values())

    # tulba asukoht X-teljel
    r1 = np.arange(len(bars1))
    r2 = [x + bar_width for x in r1]
    r3 = [x + bar_width for x in r2]
    r4 = [x + bar_width for x in r3]
    r5 = [x + bar_width for x in r4]

    # tulpade joonistamine
    ax0.bar(r1, bars1, color='c', width=bar_width, edgecolor='white', label='Geeniekspressioon')
    ax0.bar(r2, bars2, color='salmon', width=bar_width, edgecolor='white', label='RNA splaissimine')
    ax0.bar(r3, bars3, color='orchid', width=bar_width, edgecolor='white', label='Promootor')
    ax0.bar(r4, bars4, color='lightgreen', width=bar_width, edgecolor='white', label='Terminaator')
    ax0.bar(r5, bars5, color='gold', width=bar_width, edgecolor='white', label='Transkript')

    ax0.set_title('Geenide arv ühendatud komponendis vastavalt uuritud protsessile, cs <=' + str(cs_suurus))

    # telgede pealkirjadA
    ax0.set_xlabel('Geenide arv klikis')
    ax0.set_ylabel('Klikkide arv')

    # grupeeritud tulpade nimetused
    ax0.set_xticks([r + (bar_width / 5) for r in range(len(bars1))])
    ax0.set_xticklabels(["üks", "mitu"])

    ax0.legend()
    plt.show()

joonista(ge_ülekatteid, txrev_ülekatteid, upstream_ülekatteid, downstream_ülekatteid, tx_ülekatteid)
