import pandas as pd
import pyranges as pr
import numpy as np
import matplotlib.pyplot as plt
import collections

pd.set_option('display.max_columns', None)
desired_width = 1000
pd.set_option('display.width', desired_width)
pd.set_option('display.max_colwidth', 100)
np.set_printoptions(linewidth=desired_width)

'''
Failist andmete sisselugemine, veergude lisamine ja kustutamine.
Sisend: failinimi
Väljund: tabel sobilike veergudega Pyranges ja muu töötluse jaoks
'''
def algtöötlus(failinimi):
    tabel = pd.read_csv(r"C:\Users\pr\Documents\kool\6 semester 2020-21\Baka\andmed\\" + failinimi, sep="	")  # loeb faili sisse
    tabel = tabel[['phenotype_id', 'chr', 'cs_id', 'pip', 'cs_size', 'variant_id']]  # mis veerud jätta alles
    tabel = tabel.rename(columns={'chr': 'Chromosome'})
    tabel = pd.concat([tabel, tabel["variant_id"].str.split("_", expand=True)[1].rename("Start")], axis=1)  # uus veerg Start
    tabel["Start"] = pd.to_numeric(tabel["Start"])
    tabel['End'] = tabel['Start'] + 1  # Pyranges puhul peab mingi vahemik olema, suurema vahemiku puhul loetakse liiga palju ülekatteid
    tabel = tabel[['phenotype_id', 'cs_id', 'Chromosome', 'Start', 'End', 'pip', 'cs_size']]  # mis veerud jätta alles

    return tabel


'''
Tabelist sobilike credible set-ide välja valimine
Sisend: tabel, credible set-i maksimaalne suurus, kas on txrev või ge fail
Väljund: pyranges tabel, kus on iga geeni kohta 1 sobilik cs
'''
def cs_filtreerimine(tabel, cs_suurus, onSplaissimine):
    tabel_cs = tabel.copy()
    tabel_cs = tabel_cs[tabel_cs['cs_size'] < cs_suurus]  # credible set maksimaalne suurus

    if (onSplaissimine):  # vaid txrev puhul
        tabel_cs = tabel_cs[tabel_cs['phenotype_id'].str.contains("contained")]  # jätab alles vaid andmed splaissimise kohta
        tabel_cs['phenotype_id'] = tabel_cs['phenotype_id'].str.split(".", expand=True)  # jätab alles vaid geeni nime

    tabel_cs = tabel_cs.sort_values('pip', ascending=False).drop_duplicates(
        ['phenotype_id'])  # jätab iga geeni kohta vaid ühe CS-i alles, mille pip väärtus on kõige suurem

    # sobilikest cs-idest teeb Pyranges tabeli
    tabel_dict = tabel_cs.to_dict()
    cs_pr_tabel = pr.from_dict(tabel_dict)

    return cs_pr_tabel


'''
Tabelist ülekatete leidmine kasutades Pyranges count overlaps meetodit
Sisend: tabel, cs_list, faili nime algus
Väljund: leitud ülekatted
'''
def ülekatete_leidmine(tabel, cs_pr_tabel, nimi):
    tabel_dict = tabel.to_dict()
    pr_tabel = pr.from_dict(tabel_dict)  # teeb Pyranges tabeli

    # loendab kogu tabeli ülekatteid varem välja filtreeritud sobilikest geenidest, salvestab tulemuse NumberOverlaps veergu
    overlaps = pr_tabel.count_overlaps(cs_pr_tabel)
    overlaps = overlaps[overlaps.NumberOverlaps > 0]  # jätab alles vaid vähemalt ühe ülekatte leidnud cs-i
    ülekatteid = overlaps.NumberOverlaps.value_counts().to_dict()  # salvestab leitud ülekatete arvud sõnastikku

    # overlaps.rp()
    # print(overlaps)
    overlaps.to_csv(path=nimi + "_overlaps.csv")  # salvestab leitud vähemalt ühe ülekattega variandid faili

    return ülekatteid


'''
Ülekatete sõnastike sama pikkuseks tegemine ja järjestamine, et oleks lihtsam joonist teha
Sisend: ge ja txrev ülekatete sõnastikud
Väljund: korrastatud ge ja txrev sõnastikud
'''
def ülekatete_korrastamine(ge_count, txrev_count):
    suurim_pikkus = max(len(ge_count), len(txrev_count))

    for sõnastik in (ge_count, txrev_count):
        while (len(sõnastik) < suurim_pikkus):
            sõnastik[len(sõnastik) + 1] = 0

    ge_count = collections.OrderedDict(sorted(ge_count.items()))  # järjestab sõnastiku
    txrev_count = collections.OrderedDict(sorted(txrev_count.items()))  # järjestab sõnastiku

    return ge_count, txrev_count


cs_suurus = 10  # credible set suurus (kasutatakse filtreerimises), mis peab sellest jääma väiksemaks

# geeniekspressioon
print("\nGeeniekspressioon:")
# ge_tabel = algtöötlus("ge.txt")
ge_tabel = algtöötlus("HipSci.iPSC_ge.purity_filtered.txt")
ge_cs_list = cs_filtreerimine(ge_tabel, cs_suurus, False)
ge_count = ülekatete_leidmine(ge_tabel, ge_cs_list, 'ge')
print("Geeniekspressiooniga leitud ülekatted: ", ge_count)

# RNA splaissimine
print("\nRNA splaissimine:")
txrev_tabel = algtöötlus("HipSci.iPSC_txrev.purity_filtered.txt")
# txrev_tabel = algtöötlus("txrev.txt")
txrev_cs_list = cs_filtreerimine(txrev_tabel, cs_suurus, True)
txrev_count = ülekatete_leidmine(txrev_tabel, txrev_cs_list, 'txrev')
print("RNA splaissimisega leitud ülekatted: ", txrev_count)

ge_count, txrev_count = ülekatete_korrastamine(ge_count, txrev_count)


def joonista(ge, txrev):
    fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(16, 8))

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
    ax0.set_ylabel('Geenide arv')

    ax0.set_xticks([r + (barWidth / 2) for r in range(len(bars1))])
    ax0.set_xticklabels(list(ge_count.keys()))
    for p in ax0.patches:
        ax0.annotate(str(p.get_height()), (p.get_x() * 1.005 + barWidth / 6, p.get_height() * 1.005))

    ax0.legend()

    ge.pop(1, None)
    txrev.pop(1, None)

    # set height of bar
    bars1 = list(ge.values())
    bars2 = list(txrev.values())

    # Set position of bar on X axis
    r1 = np.arange(len(bars1))
    r2 = [x + barWidth for x in r1]

    # Make the plot
    ax1.bar(r1, bars1, color='c', width=barWidth, edgecolor='white', label='Geeniekspressioon')
    ax1.bar(r2, bars2, color='coral', width=barWidth, edgecolor='white', label='RNA splaissimine')

    ax1.set_title('Tulemuste võrdlus 2+ ülekattega')

    # Add xticks on the middle of the group bars
    ax1.set_xlabel('Ülekatteid')
    ax1.set_ylabel('Geenide arv')
    ax1.set_xticks([r + (barWidth / 2) for r in range(len(bars1))])
    ax1.set_xticklabels(list(ge_count.keys()))

    for p in ax1.patches:
        ax1.annotate(str(p.get_height()), (p.get_x() * 1.005 + barWidth / 6, p.get_height() * 1.005))

    # Create legend & Show graphic
    ax1.legend()
    plt.show()


joonista(ge_count, txrev_count)
