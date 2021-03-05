import pandas as pd
import numpy as np
import pyranges as pr
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt


pd.set_option('display.max_columns', None)
desired_width = 1000
pd.set_option('display.width', desired_width)
pd.set_option('display.max_colwidth', 100)
np.set_printoptions(linewidth=desired_width)


'''
Failist andmete sisselugemine, veergude lisamine ja kustutamine.
Sisend: failinimi
Väljund: tabel sobilike veergudega
'''
def algtöötlus(failinimi):
    tabel = pd.read_csv(r"C:\Users\pr\Documents\kool\6 semester 2020-21\Baka\andmed\\" + failinimi, sep="	")
    tabel['full_pos'] = tabel['chr'].astype(str) + "-" + tabel['pos'].astype(str) # uus veerg koos kromosoomi ja asukohaga
    tabel = tabel[['phenotype_id', 'chr', 'full_pos', 'cs_id', 'pip', 'cs_size']]  # mis veerud jätta alles
    return tabel


'''
Tabelist sobilike credible set-ide välja valimine
Sisend: tabel, credible set-i maksimaalne suurus, kas on txrev või ge fail
Väljund: list credible set-idest, mis vastavad tingimustele
'''
def cs_filtreerimine(tabel, cs_suurus, onSplaissimine):
    tabel_cs = tabel.copy()
    tabel_cs = tabel_cs[tabel_cs['cs_size'] < cs_suurus]  # credible set maksimaalne suurus

    if (onSplaissimine): # vaid txrev puhul
        tabel_cs = tabel_cs[tabel_cs['phenotype_id'].str.contains("contained")]  # jätab alles vaid andmed splaissimise kohta
        tabel_cs['phenotype_id'] = tabel_cs['phenotype_id'].str.split(".", expand=True)  # jätab alles vaid geeni nime

    # TODO: kas jätta nii, et järgmine rida jätab alles vaid ühe cs-i iga geeni kohta
    tabel_cs = tabel_cs.sort_values('pip', ascending=False).drop_duplicates(['phenotype_id'])  # jätab iga geeni kohta vaid ühe CS-i alles, mille pip väärtus on kõige suurem

    cs_list = tabel_cs['cs_id'].tolist()  # sobilikud CS-id salvestab listi

    #print("Enne duplikaatide eemaldamist", len(cs_list))
    #cs_list = list(dict.fromkeys(cs_list)) #eemaldab duplikaadid
    #print("Pärast duplikaatide eemaldamist", len(cs_list))
    return cs_list


'''
Tabelist ülekatete leidmine sõnastiku abil
Sisend: tabel, cs_list
Väljund: leitud ülekatted, 
'''
def ülekatete_leidmine(tabel, cs_list):
    print("Algne tabeli suurus:                    ", len(tabel))
    tabel = tabel[tabel.cs_id.isin(cs_list)].copy()  # jätab alles vaid sellised veerud, milles cs_id on juba varem välja valitud
    print("Tabeli suurus ainult sobilike cs-idega: ", len(tabel))

    cs_pos_list = tabel['full_pos'].tolist()
    asukoht_cs_sõnastik = defaultdict(list)

    for i in range(len(tabel)):
        if (tabel.iloc[i]['full_pos'] in cs_pos_list):
            cs = tabel.iloc[i]['cs_id']
            asukoht = tabel.iloc[i]['full_pos']
            asukoht_cs_sõnastik[asukoht].append(cs) #TODO: mille põhjal eristada CS-i?, hetkel on kogu 'cs-id' informatsioon oluline (seega sama geeni kohta

    #print("dictionary_list", asukoht_cs_sõnastik)

    ülekatteid_kokku = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0}
    #ülekatteid_kokku = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0, 12: 0, 13: 0, 14: 0} #mingitel tingimustel tuleb palju ülekatteid
    #ülekatteid_kokku = defaultdict(int)

    for i in asukoht_cs_sõnastik:
        ülekatteid = len(asukoht_cs_sõnastik.get(i))
        ülekatteid_kokku[ülekatteid] += 1

    print("Ülekatted: ", ülekatteid_kokku)

    kokku_ridu = len(asukoht_cs_sõnastik)
    return ülekatteid_kokku, kokku_ridu

cs_suurus = 10
# geeniekspressioon
print("\nGeeniekspressioon:")
ge_tabel = algtöötlus("HipSci.iPSC_ge.purity_filtered.txt")
#ge_tabel = algtöötlus("ge.txt")
ge_cs_list = cs_filtreerimine(ge_tabel, cs_suurus, False)
ge_count, ge_summa = ülekatete_leidmine(ge_tabel, ge_cs_list)

print("\nRNA splaissimine:")
txrev_tabel = algtöötlus("HipSci.iPSC_txrev.purity_filtered.txt")
#txrev_tabel = algtöötlus("txrev.txt")
txrev_cs_list = cs_filtreerimine(txrev_tabel, cs_suurus, True)
txrev_count, txrev_summa = ülekatete_leidmine(txrev_tabel, txrev_cs_list)

print("\nGeeniekspressiooniga leitud ülekatete arv: ", ge_summa)
print("RNA splaissimisega leitud ülekatete arv: ", txrev_summa)


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
    ax1.set_xlabel('Ülekatted')
    ax1.set_ylabel('Geenide arv')
    ax1.set_xticks([r + (barWidth / 2) for r in range(len(bars1))])
    ax1.set_xticklabels(list(ge_count.keys()))

    for p in ax1.patches:
        ax1.annotate(str(p.get_height()), (p.get_x() * 1.005 + barWidth / 6, p.get_height() * 1.005))

    # Create legend & Show graphic
    ax1.legend()
    plt.show()


joonista(ge_count, txrev_count)
