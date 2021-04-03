import scipy.stats as stats

#Nullhüpotees - RNA splaissimise ja geeniekspressiooni puhul on ühe geeni ülekatete suhe mitmega sama
#kui p-väärtus on väiksem kui 0,05, siis saame nullhüpoteesi kummutada

oddsratio, pvalue = stats.fisher_exact([[1878, 81], [1863, 16]])
print(oddsratio) # hinnang 95%-usaldusintervall šansside suhetele
print(pvalue) #kas on statistiliselt oluline (statistiline seos)

